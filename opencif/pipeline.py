"""opencif.pipeline: a specification in, and a photomask, its renders, and a validation record out.

    run(spec)  -> record (dict); writes <out>/<name>.dxf, .json, .pkl, _summary.md, _layout/_zooms/_physics.png

A run validates the spec, applies the seal policy for meander heads, evaluates the design gates in every mode, fits the
layout to the chip area, builds and measures the drawn mask, assesses every check, and renders the figures. In more
detail, run() first clears the engine's caches, loads and validates the spec, and removes the files that an earlier run
of the same name left in the output folder, so that a record never mixes runs. It then plans the device (prepare),
refuses it before anything is drawn when a refusal applies (refusal), fits the layout to the chip (fit_layout), draws
and validates the mask (build.build), assesses every check (assess), renders the figures, names the photomask after its
verdict (name_mask), composes the array when the spec asks for one (compose.array), and writes the record and its
summary.

There are three gate modes. In strict mode, a design outside the validity envelope (the kappa plateau u >= 0.035, gap
dominance M13 < 0.05, channel Re <= 10, and the Stokes gates Re_gap and Re_delta <= 1e-2) is refused, so the record
reports the gates and the cutoff band and no mask is drawn. This follows App. B's thesis that the design code refuses a
design that violates the gates (CIF_v2's M13 and M14). Strict mode also refuses a meander head whose fold width over the
gap, GS/G, lies outside gamma_GS [1.4, 1.5] (spec v5 sec. 14.4 and the sec. 18 pre-flight; CIF_v2 design()). A design
over its chip area is drawn and reported, and it is never refused or failed (App. B thesis; C8 CQ5). In evaluate mode
the mask is drawn and fully validated, and the record says which gate the design leaves. Reproduce mode works as
evaluate mode does, and it also allows and lists the reproduction-only options (the papers' deliberate departures).
"""
import json
import math
import os
import time

import numpy as np

from . import spec as SP
from . import build as BU
from .checks import auditor as AU
from .figures import proof as RN
from opencif.physics import schedule as S  # noqa: E402
from opencif.geometry import layout as L  # noqa: E402

RUNS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'runs')


def _rows_txt(rows, max_parts=6):
    """Return a compact list of rows such as '3-5, 9, 12-14', with at most max_parts runs followed by ' ...' (with
    max_parts None, every run is listed)."""
    if not rows:
        return '-'
    out, a, b = [], rows[0], rows[0]
    for r in rows[1:]:
        if r == b + 1:
            b = r
        else:
            out.append('%d' % a if a == b else '%d-%d' % (a, b))
            a = b = r
    out.append('%d' % a if a == b else '%d-%d' % (a, b))
    if max_parts is None or len(out) <= max_parts:
        return ', '.join(out)
    return ', '.join(out[:max_parts]) + ' ...'


def code_stamp():
    """Return the code stamp, a sha256 over the package's source, so that every record says which code produced it.
    The hash covers every Python file under opencif/ and its data (JSON), each hashed with its path relative to the
    package and with '/' as the separator on every system."""
    import hashlib
    import glob
    root = os.path.dirname(os.path.abspath(__file__))
    h = hashlib.sha256()
    files = [f for f in glob.glob(os.path.join(root, '**', '*.py'), recursive=True) +
             glob.glob(os.path.join(root, 'data', '*.json')) + glob.glob(os.path.join(root, 'data', '*', '*.json'))
             if '__pycache__' not in f]
    posix = lambda f: os.path.relpath(f, root).replace(os.sep, '/')
    files.sort(key=posix)                                # the same order on every system
    for f in files:
        h.update(posix(f).encode())
        h.update(open(f, 'rb').read())
    from . import __version__
    return dict(version=__version__, sha256=h.hexdigest()[:16], files=len(files))


CODE = code_stamp()


def environment_stamp():
    """Return what the code stamp does not cover, which is the test suite's hash, the lock file's hash and its pins
    against the installed versions, and the interpreter and platform (ADJ-C7). The stamp is recorded beside the code
    stamp (record 'environment') rather than inside it, so the code stamp keeps the definition it had in RC6."""
    import hashlib
    import glob
    import platform
    from importlib import metadata
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    posix = lambda f: os.path.relpath(f, root).replace(os.sep, '/')
    tests = sorted(glob.glob(os.path.join(root, 'tests', '*.py')), key=posix)
    th = hashlib.sha256()
    for f in tests:
        th.update(posix(f).encode())
        with open(f, 'rb') as fh_:
            th.update(fh_.read())
    lock_p = os.path.join(root, 'requirements-tested.txt')
    pins, lock_sha = {}, None
    if os.path.exists(lock_p):
        with open(lock_p, 'rb') as fh_:
            raw = fh_.read()
        lock_sha = hashlib.sha256(raw).hexdigest()[:16]
        # A lock file that is not UTF-8 cannot break the import (R16).
        for ln in raw.decode('utf-8', errors='replace').splitlines():
            ln = ln.split('#')[0].strip()
            if '==' in ln:
                k_, v_ = ln.split('==', 1)
                pins[k_.strip().lower()] = v_.strip()
    inst = {}
    for k_ in sorted(set(pins) | {'numpy', 'scipy', 'shapely', 'matplotlib', 'pyyaml', 'ezdxf', 'pillow'}):
        try:
            inst[k_] = metadata.version(k_)
        except metadata.PackageNotFoundError:
            inst[k_] = None
    return dict(tests_sha256=(th.hexdigest()[:16] if tests else None), tests_files=len(tests),
                lock_file='requirements-tested.txt', lock_sha256=lock_sha,
                installed=inst, matches_lock=(bool(pins) and all(inst.get(k_) == v_ for k_, v_ in pins.items())),
                python=platform.python_version(), platform=platform.platform(),
                note='rc7 (ADJ-C7): the test suite and the environment, beside the code stamp (not hashed into it)')


ENV = environment_stamp()


def is_sealed(sc):
    """Return True when the plan seals the inlet, that is, when gap 1 extracts f* (to a relative 1e-9)."""
    return abs(sc.f_seed / sc.info['fstar'] - 1.0) < 1e-9


def _plan(spec):
    """Plan the device and return (schedule, sealed, error). A plan that raises RuntimeError gives (None, False, its
    message)."""
    try:
        sc = S.plan(spec)
        return sc, is_sealed(sc), None
    except RuntimeError as e:
        return None, False, str(e)


def seal_policy(spec, opts, notes, say):
    """Apply the seal policy of a meander head and return the spec.

    A sealed inlet (gap 1 at f*, so that its cutoff is D) needs room for the head's long first units. In the matrix
    design the first unit is 24 mm of 28 um channel, and each later unit is about 1/k of it. The head style asked for is
    the head drawn, and no other head is ever substituted. The three policies are:

      'require'  sealed, or the build is refused with the numbers (the arms needed against the room for them, the head
                 width that would seal, and whether the narrow head would)
      'auto'     sealed if the asked head can seal; if it cannot, gap 1 is allowed to extract more outside strict mode
                 (declared and reported, as under 'allow'), and in strict mode it is refused as under 'require'
      'allow'    gap 1 may extract more than f* (its cutoff is reported and the row declared); not in strict mode

    A refusal is left in opts['seal_refusal'] for `run` to record."""
    import copy
    sc, ok, err = _plan(spec)
    if ok:
        return spec
    style = opts['style']
    if sc is None:                                       # not even an unsealed head fits; refused with the reason
        opts['seal_refusal'] = ('seal: the %s head cannot be planned within the %.0f um head width (%s); a wider head '
                                'or chip, or another style, is needed' % (style, spec.W, err))
        say(opts['seal_refusal'])
        return spec
    D_um = S.KAPPA_G * spec.delta_star
    raised = 'gap 1 extracts %.2f f* (cutoff %.3f um against D %.3f um)' % (
        sc.f_seed / sc.info['fstar'], S.KAPPA_G * S.x_of(sc.f_seed) * spec.wc0, D_um)
    if opts['seal'] == 'allow' or (opts['seal'] == 'auto' and opts['mode'] != 'strict'):
        notes.append('seal: the %s head does not seal gap 1 within the %.0f um head width: %s; allowed by meander.seal: '
                     '%s (gates.mode %s), and gap 1 is a declared row' % (style, spec.W, raised, opts['seal'], opts['mode']))
        say(notes[-1])
        opts['seal'] = 'allow'                            # the seed row is declared like an explicit allow
        return spec
    need, cap = sc.info.get('A_need_sealed'), sc.info.get('A_cap')
    W_seal = None
    for Wt in np.linspace(spec.W * 1.1, spec.W * 3.0, 12):
        spt = copy.copy(spec)
        spt.W = float(Wt)
        if _plan(spt)[1]:
            W_seal = float(Wt)
            break
    narrow_ok = False
    if style != 'narrow':
        spn = copy.copy(spec)
        spn.head, spn.w_narrow, spn.w_arm = 'pack', spec.F_min, spec.F_min
        narrow_ok = _plan(spn)[1]
    msg = ('seal: the %s head cannot seal gap 1 within the %.0f um head width%s: unsealed, %s.  Remedies: %s%s; or '
           'meander.seal: allow in gates.mode: evaluate (gap 1 declared).' % (
               style, spec.W, (' (it needs %s arms; %s fit)' % (need, cap)) if (need is not None and cap is not None
                                                                                and style in ('pack', 'arms', 'narrow')) else '',
               raised,
               ('a head about %.0f um wide (a chip at least %.1f mm wide) seals it' % (
                   W_seal, math.ceil((W_seal + 2 * opts['margin']) / 100.0) / 10.0))
               if W_seal else 'no head up to 3 times as wide seals it',
               '; the narrow head (process-floor units) would seal it at this width' if narrow_ok else ''))
    opts['seal_refusal'] = msg
    say(msg)
    return spec


def shrink_head(spec, opts, notes, say):
    """Return the spec with the narrowest head width W that still seals gap 1 (+3 %), searched from 0.3 to 1 times the
    available width in every layout. The packing head then fills a compact head the way Dinh's does (arms reaching
    upstream, a few lane pairs each), rather than folding its first arm into a tall comb of short lanes across the
    whole chip width (m01 on a 12 mm strip, 37 lane pairs in 10 mm). In folded layouts it also keeps the leg pitch
    small. A head width given in the spec is kept."""
    import copy
    if not opts['head_width_auto'] or spec.head == 'none':
        return spec
    if not _plan(spec)[1]:
        return spec
    lo, hi = 0.3 * spec.W, spec.W
    spt = copy.copy(spec)
    spt.W = lo
    if _plan(spt)[1]:
        hi = lo
    else:
        for _ in range(7):
            mid = 0.5 * (lo + hi)
            spt = copy.copy(spec)
            spt.W = mid
            if _plan(spt)[1]:
                hi = mid
            else:
                lo = mid
    W_new = math.ceil(hi * 1.03 / 50.0) * 50.0
    spt = copy.copy(spec)
    spt.W = W_new
    if not _plan(spt)[1]:                                # sealing need not be monotone in W, so the given width is kept
        notes.append('head width: kept %.0f um (the narrowest-width search did not give a width that seals)' % spec.W)
        return spec
    if W_new < spec.W:
        notes.append('head width: %.0f um, the narrowest that seals gap 1 within 0.3 to 1 times the available %.0f um '
                     '(+3%%)%s' % (W_new, spec.W, '; 0.3 times seals already: the floor of the search' if hi == lo else ''))
        say(notes[-1])
        spec = copy.copy(spec)
        spec.W = W_new
    return spec


def floor_start(spec, opts, notes, say):
    """Apply Gifford 2014's start rule and return the spec, with the side flow S0 that the rule needs. Gifford writes
    (p. 4498) that w_s(1) 'may be made close to or equal to zero, depending on the device production capabilities
    available', so the first side channel is drawn at the stated width (by default the process floor), and the closure
    designs every later row from it. In the design this is a side flow S0 already present at gap 1. The drawn device is
    sealed, however, so its first gaps extract more than f* until the side flow catches up. The record's check of the
    realized cutoff measures that entrance error on the drawn mask."""
    import copy
    w1 = opts['start_width']
    sc0 = S.plan(spec)
    if sc0.ws[0] >= w1:
        notes.append('start: the sealed design already starts at %.2f um >= %.2f um; no floor start needed' % (sc0.ws[0], w1))
        return spec
    lo, hi = 0.0, 1e-3
    while True:
        sp = copy.copy(spec); sp.S0 = hi
        if S.plan(sp).ws[0] >= w1:
            break
        lo, hi = hi, hi * 4
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        sp = copy.copy(spec); sp.S0 = mid
        if S.plan(sp).ws[0] >= w1:
            hi = mid
        else:
            lo = mid
    sp = copy.copy(spec); sp.S0 = hi
    notes.append('start: floor rule (Gifford 2014 p. 4498): w_s(1) = %.2f um (sealed design: %.2f um); the design assumes a '
                 'side flow S0 = %.3e of the inlet flow at gap 1, which the sealed drawing does not have' % (w1, sc0.ws[0], hi))
    if opts.get('mode') != 'reproduce':
        # The entrance rows that the floor start over-extracts are declared in reproduce mode only
        # (build.active_options), so a new design with a floor start fails Q2 on them, as m14 shows. The note says so
        # before drawing (review8's stress test).
        notes.append('start: outside gates.mode: reproduce the floor start\'s entrance rows are not declared: where they '
                     'over-extract beyond checks.cutoff_tolerance (the usual case; the larger S0, the more rows) the build '
                     'FAILS Q2 on them (its photomask named <name>_FAILED.dxf); a meander head starts sealed at the floor '
                     'instead')
    say(notes[-1])
    return sp


def leg_candidates(sc, kind, opts):
    """Return the leg counts that the layout tries, in order. These are the fixed count, the default without a chip,
    or (auto, with a chip) the fewest whose unfolded length could fit and the next nine (odd only for a serpentine)."""
    if kind == 'linear':
        return [1]
    legs = opts['legs']
    if isinstance(legs, int):
        return [legs]
    chip_w, chip_l = opts['chip']
    if not (chip_w and chip_l):
        return [5 if kind == 'serpentine' else 4]
    Lc = chip_l * 1000.0 - 2 * opts['margin']
    path = float(sc.x[-1]) + 3000.0                       # x runs from gap 1, so the unfolded head is included
    n0 = max(2 if kind == 'paperclip' else 3, int(math.ceil(path / Lc)))
    if kind == 'serpentine' and n0 % 2 == 0:
        n0 += 1
    step = 2 if kind == 'serpentine' else 1
    return list(range(n0, n0 + 10 * step, step))


def outlet_precheck(spec, sc, opts, notes, say=lambda m: None):
    """Check whether the requested outlet can close the design at all.

    The outlet closes the design by its lead's resistance, which is planned with the layout's last leg, because the
    leg's straight stretch to the outlet moves with the leg count. This function plans the layout, without drawing, for
    the leg counts that the fitter would try (`leg_candidates`). If the outlet cannot close at any of them, the design
    is refused before drawing in every mode (opts['outlet_refusal'], exit 5), as `check` and generate(check_only=True)
    also report. In strict mode an outlet whose stretch-end term lies beyond the solved 2-D tables is refused here as
    well (layout.OutletOutsideTables), whereas the other modes build it and say so in the verdict. A count that fails
    for another reason (too short for n legs) ends the search, as it ends the fitter's. The fitter repeats the outlet
    plan for a widened head."""
    kind = opts['layout']
    errs, first = [], None
    for i, n in enumerate(leg_candidates(sc, kind, opts)):
        try:
            L.plan_layout(sc, kind, n, opts.get('leg_spacing'))
            return
        except L.OutletOutsideTables as e:                # the stretch-end tables' domain (strict mode)
            opts['outlet_refusal'] = ('the %s outlet: %s: build it in evaluate or reproduce mode, where the end values '
                                      'are held and reported, or change the outlet' % (opts.get('outlet'), e))
            say(opts['outlet_refusal'])
            return
        except L.OutletInfeasible as e:
            errs.append('%d leg(s): %s' % (n, e))
        except (ValueError, RuntimeError) as e:
            if i == 0:                                    # the fitter's first count cannot be laid out, so none can
                first = '%d leg(s): %s' % (n, e)
            break
    if first and not errs:                                # a spec error before drawing, as the fitter would give
        opts['layout_problem'] = ('the %s layout cannot be laid out at the leg count%s the fitter tries (%s); more legs only '
                                  'shorten each, so none can: use a linear layout or a longer device'
                                  % (kind, '' if isinstance(opts['legs'], int) else 's', first))
    if errs:
        opts['outlet_refusal'] = ('the %s outlet cannot close this design (F %.4f) at any leg count the layout tries: %s'
                                  % (opts.get('outlet'), spec.F_target, errs[0] if len(errs) == 1 else
                                     '%s (and %d more counts alike)' % (errs[0], len(errs) - 1)))
        say(opts['outlet_refusal'])


def _fit_once(spec, sc, kind, opts, Lc, Wc, say):
    """Find the fewest legs that fit (Lc x Wc inside the margin) by drawing each candidate, and return
    (n or None, best, tried)."""
    path = float(sc.x[-1]) + 3000.0                       # x runs from gap 1, so the unfolded head is included
    if sc.NM > 0:
        # The head is drawn unfolded at the start of the first leg. If the head alone is longer than the chip, no leg
        # count can fit, so the futile drawings are skipped and the caller may widen the head.
        head_len = float(sc.x[sc.NM]) - float(L.head_top_x(sc))
        if head_len > Lc:
            say('fit: the %.0f um head alone is %.2f mm long, more than the chip allows (%.2f mm)' % (spec.W, head_len / 1e3, Lc / 1e3))
            return None, None, ['head %.0f um: %.2f mm long alone' % (spec.W, head_len / 1e3)]
    n0 = max(2 if kind == 'paperclip' else 3, int(math.ceil(path / Lc)))
    if kind == 'serpentine' and n0 % 2 == 0:
        n0 += 1
    step = 2 if kind == 'serpentine' else 1
    best, tried, collided, curable = None, [], None, None
    for n in range(n0, n0 + 10 * step, step):
        try:
            sc_l, legs_ = L.plan_layout(sc, kind, n, opts.get('leg_spacing'))
            fl, sc_d, info_d = L.assemble(sc_l, L.make_legs(sc_l))
        except L.OutletInfeasible as e:
            # The outlet cannot close at n legs, because its stretch moves with n, so the next count is tried.
            # outlet_precheck refuses the design when no count can close it.
            tried.append('%d legs: %s' % (n, e))
            say('fit: %d legs: %s; next' % (n, e))
            continue
        except (ValueError, RuntimeError) as e:          # too short for n legs, or loops past the end; n stops growing
            tried.append('%d legs: %s' % (n, e))
            say('fit: %d legs not drawable (%s)' % (n, e))
            break
        # A candidate whose parts collide is rejected. A two-leg paperclip ends beside its inlet, and in m10 the
        # outlet's filtrate port overlapped the inlet port. A candidate must therefore be one polygon with the planned
        # holes and separated ports, with every part at least the floor from every part it is not joined to, and with
        # the outlet ports clear of the inlet port.
        single, holes, exp_holes = BU.topology_counts(fl, sc_d, info_d)
        sep = BU.port_separation(fl, sc_d, info_d)['separated'] if single else False
        pc = L.part_clearance(info_d['parts'], spec.F_min)
        ps = L.port_spacing(info_d, spec)
        if not (single and len(holes) == exp_holes and sep and pc['ok'] and ps['ok']):
            if curable is None and kind != 'linear' and (not ps['ok'] or not sep or not pc['ok']):
                # A count whose only collision is its outlet meeting the inlet end is one that the build cures by
                # lengthening the outlet stretch (build.plan_and_draw). When every count collides, such a count is
                # preferred over the first count that collided (code review 4, F7).
                try:
                    if BU._outlet_meets_inlet(BU.collisions(fl, sc_d, info_d, spec)):
                        curable = n
                except Exception:                         # a drawing that the checks cannot read is not curable
                    pass
            why = ('%d fluid pieces' % len(fl.geoms) if not single else
                   'holes %d of %d' % (len(holes), exp_holes) if len(holes) != exp_holes else
                   'outlet ports not separated' if not sep else
                   '%s and %s %.1f um apart' % (pc['closest'][0]['parts'][0], pc['closest'][0]['parts'][1],
                                                pc['closest'][0]['gap_um']) if not pc['ok'] else
                   'outlet port %s %.0f um from the inlet port (at least %.0f)' % (
                       ps['port'], ps['inlet_to_nearest_outlet_port_um'], ps['required_um']))
            tried.append('%d legs: parts collide (%s)' % (n, why))
            say('fit: %d legs: parts collide (%s), skipped' % (n, why))
            b = fl.bounds
            if collided is None or curable == n:
                collided = (n, b[2] - b[0], b[3] - b[1], True)       # the fourth field marks a colliding drawing
            if b[3] - b[1] > Wc:                          # more legs only widen it further
                break
            continue
        b = fl.bounds
        lx, ly = b[2] - b[0], b[3] - b[1]
        tried.append('%d legs: %.2f x %.2f mm' % (n, lx / 1e3, ly / 1e3))
        say('fit: %d legs, head width %.0f um -> %.2f x %.2f mm (chip %.2f x %.2f mm inside the margin)'
            % (n, spec.W, lx / 1e3, ly / 1e3, Lc / 1e3, Wc / 1e3))
        if lx <= Lc and ly <= Wc:
            return n, (n, lx, ly, False), tried
        if best is None or max(lx / Lc, ly / Wc) < max(best[1] / Lc, best[2] / Wc):
            best = (n, lx, ly, False)
        if ly > Wc:                                       # more legs only widen it further
            break
    # No count fits. When some drawing does not collide, the one with the smallest overshoot of the chip is returned.
    # When every drawable count collided, the first count whose outlet only meets the inlet end is returned (the build
    # lengthens its outlet stretch), or else the first count. If that drawing still collides, the build refuses it on
    # its outline before anything is written (build.collisions) and names the collision and the remedies.
    return None, (best if best is not None else collided), tried


def fit_layout(spec, sc, opts, notes, say):
    """Choose the number of legs, which is either the number given or ('auto') the fewest that fit the chip (length
    along the legs, width across), and return (kind, legs, spec, sc).

    The device's long axis runs along the chip's length. A leg count that cannot be drawn, or whose drawn parts collide
    (the drawing must be one polygon with the planned holes, with the outlet ports separated from each other and from
    the inlet), is skipped, and a note says so. A meander head is drawn unfolded at the start of the first leg, so its
    length alone can exceed the chip. When no leg count fits and the head width is automatic, a packed head is widened
    up to the chip width, which folds its arms into more lanes and gives a shorter head and a wider first leg, and the
    fit is repeated. The returned spec and sc change only when the head was widened."""
    import copy
    kind, legs = opts['layout'], opts['legs']
    if kind == 'linear':
        return kind, 1, spec, sc
    chip_w, chip_l = opts['chip']
    if isinstance(legs, int):
        return kind, legs, spec, sc
    if not (chip_w and chip_l):
        legs = 5 if kind == 'serpentine' else 4
        notes.append('layout.legs auto without a chip area: %d legs (give chip.width_mm and chip.length_mm to fit)' % legs)
        return kind, legs, spec, sc
    Lc, Wc = chip_l * 1000.0 - 2 * opts['margin'], chip_w * 1000.0 - 2 * opts['margin']
    n, best, tried = _fit_once(spec, sc, kind, opts, Lc, Wc, say)
    if n is not None:
        notes.append('layout: %d legs, the fewest that fit the chip (%.2f x %.2f mm)' % (n, best[1] / 1e3, best[2] / 1e3))
        return kind, n, spec, sc
    best_sp = (spec, sc)
    # Widening shortens only a packed head, whose arms fold into more lanes. Hairpin units and Dinh's single-lane arms
    # have lengths fixed by their resistances, so a wider head draws the same device (m05 and m12 made three futile
    # drawings).
    if opts.get('head_width_auto') and spec.head == 'pack' and sc.A > 0:
        W = spec.W
        while W < Wc:
            W = min(Wc, W * 1.5)
            sp2 = copy.copy(spec)
            sp2.W = W
            sc2, ok, err = _plan(sp2)
            if not ok:
                continue
            n2, best2, tried2 = _fit_once(sp2, sc2, kind, opts, Lc, Wc, say)
            tried += ['head %.0f um, %s' % (W, t) for t in tried2]
            if n2 is not None:
                notes.append('layout: no leg count fitted with the %.0f um head (its unfolded length alone exceeds the chip); '
                             'head widened to %.0f um: %d legs fit (%.2f x %.2f mm)' % (spec.W, W, n2, best2[1] / 1e3, best2[2] / 1e3))
                say(notes[-1])
                return kind, n2, sp2, sc2
            # A clean drawing always ranks before a colliding one, and the smaller overshoot of the chip ranks next.
            rank = lambda b_: (b_[3], max(b_[1] / Lc, b_[2] / Wc))
            if best2 is not None and (best is None or rank(best2) < rank(best)):
                best, best_sp = best2, (sp2, sc2)
    if best is None:
        if tried and all('outlet' in t for t in tried):
            raise L.OutletInfeasible('no leg count closes the %s outlet (%s)' % (opts.get('outlet'), '; '.join(tried)))
        raise SP.SpecError('layout: no drawable leg count for %s (%s)' % (kind, '; '.join(tried)))
    if best[3]:
        notes.append('layout: every drawable leg count collides (%s); %d legs (%.2f x %.2f mm) go to the build, which '
                     'lengthens the outlet stretch if only the outlet meets the inlet end, or refuses the outline before '
                     'anything is written (as named above)' % (
                         '; '.join(tried), best[0], best[1] / 1e3, best[2] / 1e3))
    else:
        notes.append('layout: NO leg count fits the chip (%s); using %d legs (%.2f x %.2f mm): the record reports it '
                     'over the chip area (a packaging budget, reported, not a pass condition: App. B C8 CQ5)'
                     % ('; '.join(tried), best[0], best[1] / 1e3, best[2] / 1e3))
    if best_sp[0] is not spec:
        notes.append('head width: %.0f um (the closest fit)' % best_sp[0].W)
    return kind, best[0], best_sp[0], best_sp[1]


def requested_vs_drawn(spec, opts, cd):
    """Compare the drawing's components (build.components_drawn) with the request and return a list of mismatches."""
    bad = []
    style = opts.get('style')
    widths = cd['head_unit_widths_um']
    if style is None:
        if cd['head_units']:
            bad.append('a meander head was drawn but none was asked for')
    else:
        if not cd['head_units']:
            bad.append('no meander head was drawn (asked: %s)' % style)
        full = round(float(spec.GS), 3)
        if style == 'hairpins' and (cd['head_arms'] or any(abs(w - full) > 1e-6 for w in widths)):
            bad.append('hairpins asked; drawn: %d arms, unit widths %s' % (cd['head_arms'], widths))
        if style == 'arms':
            n_ask = int(getattr(spec, 'n_arms_fixed', 0) or 0)
            if cd['head_arms'] != n_ask or any(n != 1 for n in cd['arm_lane_pairs']):
                bad.append("Dinh's arms asked (%d single-lane); drawn: %d arms, lane pairs %s" % (
                    n_ask, cd['head_arms'], cd['arm_lane_pairs']))
        if style == 'pack' and any(abs(w - full) > 1e-6 for w in widths):
            bad.append('the packing head asked; drawn unit widths %s (channel %.1f um)' % (widths, full))
        if style == 'narrow' and getattr(spec, 'w_narrow', None) and min(widths) > float(spec.w_narrow) + 1e-6:
            bad.append('the narrow head asked; no unit at %.1f um (widths %s)' % (spec.w_narrow, widths))
        if opts.get('seal') == 'require' and not cd['gap1_sealed']:
            bad.append('seal: require, but gap 1 extracts %.3f f*' % cd['gap1_f_over_fstar'])
    if opts.get('inlet') in ('filter', 'dinh_filter') and cd['filter_posts_counted'] < 1:
        bad.append('a filter inlet asked; no filter posts drawn')
    if opts.get('inlet') == 'plain' and cd['filter_posts_counted']:
        bad.append('a plain inlet asked; %d posts drawn in it' % cd['filter_posts_counted'])
    if opts.get('outlet') and cd['outlet_style'] != opts['outlet']:
        bad.append('outlet %s asked; %s drawn' % (opts['outlet'], cd['outlet_style']))
    if cd['outlet_ports'] is not None and cd['outlet_ports'] != (3 if opts.get('outlet') == 'separate' else 2):
        bad.append('%d outlet ports drawn for a %s outlet' % (cd['outlet_ports'], opts.get('outlet')))
    if opts.get('layout') == 'linear' and cd['legs'] != 1:
        bad.append('linear layout asked; %d legs drawn' % cd['legs'])
    if isinstance(opts.get('legs'), int) and cd['legs'] != opts['legs']:
        bad.append('%d legs asked; %d drawn' % (opts['legs'], cd['legs']))
    if cd['u_turns'] != cd['legs'] - 1:
        bad.append('%d legs but %d U-turns' % (cd['legs'], cd['u_turns']))
    return bad


def band_txt(band):
    """Return the predicted cutoff band in words (build.gates' cutoff_band_um)."""
    if not band.get('band_um'):
        return 'no cutoff band (%s)' % band.get('basis', '')
    return 'cutoff band %.3f-%.3f um, %s' % (band['band_um'][0], band['band_um'][1], band.get('basis', ''))


def design_gates(spec, sc):
    """Evaluate the design gates (build.gates) on the plan at the cutoff D = kappa_g delta*."""
    return BU.gates(spec, sc, S.KAPPA_G * spec.delta_star)


PASSAGE_NOTE = 'passage (a warning, not a check; LIMITATIONS 26)'
PASSAGE_NOTE_MAX = 'passage at the largest realised gap cutoff (a warning, not a check; LIMITATIONS 26)'


def _renote(notes, prefix, text):
    """Replace any earlier note that starts with this prefix by text (or drop it when text is empty), because a
    reassess passes the record's notes back in (code review R4)."""
    notes[:] = [n for n in notes if not n.startswith(prefix)]
    if text:
        notes.append(text)


def passage_report(widths, D, notes, D_max=None, G=None, rows_no_cutoff=0):
    """Compare every measured channel class that a particle passed by the gaps can enter with D and 2 D (ADJ-C9,
    LIMITATIONS 26).

    The classes are the gap openings, the straight side channels (their narrowest section, the start included), the
    U-turn sides, the meander head units (the transition unit among them), and the outlet leads, with the center
    channel as the reference. Rows declared or realized above D pass larger particles (code review R3). Therefore, with
    D_max (the record's largest realized gap cutoff, D_realised_max_um), the classes are also compared with the largest
    particle that a gap passes, min(D_max, the gap opening G). The rows that feed each class are not traced, so this
    comparison bounds every class. Rows extracting f >= 1/2 on a side (rows_no_cutoff, the record's
    rows_f_half_or_more) have no cutoff and pass particles up to G. The result is a report and a warning in the notes,
    and it is never a check, because opencif has no passage rule (there is no particle-wall model, and the multiple of
    D that a channel needs is not established)."""
    rows, warn, warn_m = {}, [], []
    d_pass = None
    if D_max is not None and D_max > D:
        d_pass = float(min(D_max, G)) if G else float(D_max)
    if rows_no_cutoff and G:                                  # rows at f >= 1/2 pass up to the gap (codereview2)
        d_pass = float(G)
    for k_, w_ in widths.items():
        rows[k_] = dict(width_um=round(float(w_), 3), over_D=round(float(w_) / D, 3), below_D=bool(w_ < D),
                        below_2D=bool(w_ < 2.0 * D))
        if d_pass is not None:
            rows[k_].update(over_D_pass_max=round(float(w_) / d_pass, 3), below_D_pass_max=bool(w_ < d_pass))
        if k_ != 'centre channel' and w_ < 2.0 * D:
            warn.append('%s %.2f um (%.2f D)' % (k_, w_, w_ / D))
        if d_pass is not None and k_ not in ('centre channel', 'gap opening G') and w_ < d_pass:
            warn_m.append('%s %.2f um' % (k_, w_))
    _renote(notes, PASSAGE_NOTE, 'passage (a warning, not a check; LIMITATIONS 26): narrower than 2 D = %.2f um: %s'
            % (2.0 * D, '; '.join(warn)) if warn else None)
    if rows_no_cutoff and G:                                  # the note says what set d_pass (codereview3 m1, m2)
        why = '%d rows extract f >= 1/2 (no cutoff) and pass particles up to the gap %.2f um' % (int(rows_no_cutoff), float(G))
    elif d_pass is not None:
        why = ('rows above D pass particles up to %s = %.2f um' % (
            ('min(D_realised_max %.2f um, gap %.2f um)' % (D_max, G)) if G else ('D_realised_max %.2f um' % D_max), d_pass))
    else:
        why = ''
    _renote(notes, PASSAGE_NOTE_MAX, ('%s: %s; narrower than that: %s (the rows feeding each class are not traced)'
                                      % (PASSAGE_NOTE_MAX, why, '; '.join(warn_m))) if warn_m else None)
    return dict(D_um=round(D, 4), classes=rows, below_D=[k_ for k_, r_ in rows.items() if r_['below_D']],
                below_2D=[k_ for k_, r_ in rows.items() if r_['below_2D']],
                D_realised_max_um=(round(float(D_max), 4) if D_max is not None else None),
                D_pass_max_um=(round(d_pass, 4) if d_pass is not None else None),
                below_D_pass_max=[k_ for k_, r_ in rows.items() if r_.get('below_D_pass_max')],
                rows_f_half_or_more=int(rows_no_cutoff or 0),
                note='rc7 (ADJ-C9; code review R3): the narrowest measured width of each channel class against D, 2 D and '
                     '(when rows are realised above D) the largest particle a gap passes, min(D_realised_max, G), or G '
                     'when rows extract f >= 1/2 (no cutoff); a report and a warning in the notes, not a check (no '
                     'passage rule); the inlet stem and filter are not measured (the DRC at the floor only)')


# STATUS_CODES holds the enumerated values of assessment.status_codes, with one code per assessment.status field,
# while the status itself keeps its sentences (ADJ-C6).
STATUS_CODES = dict(hard_checks=('PASS', 'FAIL', 'INCOMPLETE'), validity_envelope=('INSIDE', 'OUTSIDE'),
                    chip_area=('NOT_GIVEN', 'FITS', 'OVER'), aspect_ratio_of_CIF_channels=('WITHIN_LIMIT', 'ABOVE_LIMIT'),
                    length_budget_G9=('NOT_GIVEN', 'WITHIN', 'OVER'), drawn_mask_assessment=('RUN',),
                    physics_audit=('PASS', 'PASS_OUTSIDE_ENVELOPE', 'FAIL', 'INCOMPLETE'), fabricated_or_measured=('NO',),
                    reproduction=('DESIGN_LEVEL', 'NOT_CLAIMED'), source_trace=('NAMED', 'NOT_NAMED'))


def _audit_code(pa, drc):
    """Return the physics audit's status code, which is INCOMPLETE when the audit did not run, or when its only hard
    failure is P9 and the DRC did not run, and the audit's own status otherwise."""
    if pa is None:
        return 'INCOMPLETE'
    if pa['hard_failures'] == ['P9'] and (not drc or 'error' in drc):
        return 'INCOMPLETE'
    return pa['status']


def _audit_status(pa, drc):
    """Return the physics audit's status line (assessment.status.physics_audit)."""
    code = _audit_code(pa, drc)
    if pa is None:
        return 'INCOMPLETE: the physics audit did not run'
    p8s = next((c for c in pa['checks'] if c['id'] == 'P8'), None)
    strict_txt = (' The auditor\'s own, stricter inlet limit (P8, Re <= 1; reported): %s.' % p8s['value']['result']) if p8s else ''
    if code == 'INCOMPLETE':
        return 'INCOMPLETE: the DRC did not run, so P9 could not be evaluated' + strict_txt
    if code == 'FAIL':
        return 'FAIL: %s (validation.physics_audit)' % ', '.join(pa['hard_failures']) + strict_txt
    return ('%s: the photomask file read back passes the physics auditor (CIF_v2\'s P1-P10 ported; P0: the file is '
            'the drawing; R0: its straight rows measured again by the auditor\'s own ruler; R6b: its posts are the '
            'process\'s)%s.' % (code, ' with M14 outside the validity envelope' if code == 'PASS_OUTSIDE_ENVELOPE'
                                 else '') + strict_txt)


def assess(rep, spec, opts, fluid, notes, chip_rect):
    """Assess the drawn mask and return the record's 'assessment', which holds the verdict, the exit code, every check,
    and the status.

    The assessment answers two questions on the drawn mask, and each is a hard check in every mode. Q1 asks whether the
    mask draws its design. It compares the drawn network's f at every gap with the intended network's (the closure
    resistances plus the declared deliberate deviations), and F drawn with F intended. Q2 asks whether the design
    realizes the targets, that is, whether every gap's realized cutoff lies within D +- tol and F lies within its
    tolerance of the target. The rows that a declared option (reproduce mode, or an explicitly allowed seed) is
    predicted, before measuring, to move are exempt from the cutoff test. The fabrication checks (the DRC must run,
    topology, continuity, ports, envelope, and corners) and the validity envelope complete the assessment."""
    V = rep['validation']
    out = []

    def chk(name, ok, value, need, kind='hard'):
        out.append(dict(check=name, ok=bool(ok), value=value, requirement=need, kind=kind))
    drc = V.get('drc')
    if drc and 'error' not in drc:
        chk('DRC at the process floor', drc.get('n_lost_open', 1) == 0 and drc.get('n_lost_solid', 1) == 0 and
            drc.get('n_components', 2) == 1 and drc.get('invalid_rings', 1) == 0,
            dict(lost_open=drc.get('n_lost_open'), lost_solid=drc.get('n_lost_solid'), components=drc.get('n_components'),
                 d_um=drc.get('d_um')), 'nothing narrower than %.1f um lost (tested at %.1f um); one fluid component' % (
                spec.F_min, spec.F_min - 0.1))
    else:
        chk('DRC at the process floor', False, (drc or {}).get('error', 'not run'), 'the DRC must run', 'incomplete')
    t = V['topology']
    chk('topology', t['single_component'] and t['valid'] and t['holes_ok'] and t['row_kinds_survive'],
        dict(holes=t['holes'], expected=t['holes_expected']), 'one valid polygon; every post and unit present')
    sc_ = V.get('side_continuity')
    if sc_:
        chk('side channels continuous (straight sections, every leg, both sides)', sc_['ok'] and sc_['above_floor'],
            dict(pieces=[(r['leg'], r['side'], r['pieces']) for r in sc_['legs'] if r['pieces'] != 1] or 'all 1',
                 narrowest_um=round(sc_['narrowest_um'], 3)),
            'one connected strip per leg and side; narrowest section >= %.1f um' % spec.F_min)
    chk('outlet ports separated', V['outlet_separation']['separated'], V['outlet_separation']['n_components_after_cut'],
        'product, filtrate and inlet ports in separate components once the outlet stretch is cut')
    pc = V.get('part_clearance')
    if pc is not None:
        chk('no intersections: parts keep the floor apart', pc['ok'], pc['closest'][:2],
            'every part (inlet and head, each leg, each U-turn, outlet) at least %.1f um from every part it is not joined '
            'to' % spec.F_min)
    ps = V.get('port_spacing')
    if ps is not None:
        chk('inlet and outlet ports apart', ps['ok'], ps['inlet_to_nearest_outlet_port_um'],
            'every outlet port at least %.0f um (edge to edge) from the inlet port' % ps['required_um'])
    cd = V.get('components_drawn')
    if cd is not None:
        bad = requested_vs_drawn(spec, opts, cd)
        chk('drawn components as requested', not bad, bad or 'head %s, inlet %s, outlet %s, %d leg(s)' % (
            opts.get('style') or 'none', opts.get('inlet'), opts.get('outlet'), cd['legs']),
            'the head style, sealing, inlet filter, outlet style and leg count drawn are the ones the spec asks for '
            '(nothing substituted or left empty)')
    chk('head inside its envelope', V['envelope']['ok'], round(V['envelope']['head_halfwidth_um'], 1), '<= W/2 = %.0f um' % (spec.W / 2))
    chk('corner radius', V['rounding']['ok'], round(V['rounding']['min_local_radius_um'], 3), '>= %.2f um' % spec.rmin_ok)
    uf = V.get('features_found')
    if uf is not None:                                   # a planned feature that the drawing does not contain (E5)
        chk('every planned feature found on the drawing', uf['ok'],
            dict(n_not_found=uf['n_not_found'], first=uf['not_found'][:3]) if not uf['ok'] else 'all',
            'every measurement section (head units, gaps, side segments, centre, loops, outlet) finds fluid where '
            'the plan puts the feature')
    fid = V['drawing_fidelity']
    f_no = fid.get('f_dev_max_where_no_cutoff', fid['f_dev_over_fstar_max_where_no_cutoff'])
    q1_ok = fid['cutoff_dev_max'] <= 2.5e-3 and f_no <= 5e-3 and fid['F_abs_dev'] <= 1e-4
    chk('Q1 drawing realises its design (every gap)', q1_ok,
        dict(cutoff_dev_max=float('%.3g' % fid['cutoff_dev_max']), at_row=fid['cutoff_dev_row'],
             f_dev_where_no_cutoff=float('%.3g' % f_no), f_row=fid['f_dev_row'],
             F_drawn=round(fid['F_drawn'], 6), F_intended=round(fid['F_intended'], 6)),
        'at every gap the drawn cutoff within 0.25% of the intended one, judged on the larger of the intended cutoff '
        'and D (where a row has no cutoff: f within 0.5% of the larger of |f| and f*); F within 1e-4')
    net = V['network']
    dF = abs(net['F_as_drawn'] - net['F_target'])
    chk('Q2 filtrate fraction F', dF <= opts['F_tol'], round(net['F_as_drawn'], 6), 'F = %.5f +- %g' % (net['F_target'], opts['F_tol']))
    rc = V['realised_cutoff']
    tol = opts['cut_tol']
    chk('Q2 cutoff D at every gap (undeclared rows)', rc['undeclared_rows_outside'] == 0,
        dict(rows_outside=rc['undeclared_rows_outside'], first_rows=rc['undeclared_rows_list'][:12],
             worst_over=rc['undeclared_worst_over'], worst_under=rc['undeclared_worst_under'],
             undefined_rows=rc['undeclared_undefined_rows'],
             D_range_um=[round(rc['D_realised_min_um'], 4), round(rc['D_realised_max_um'], 4)]),
        'every gap not covered by a declared option realises D +- %g%% on both sides (no reverse flow, f < 1/2)' % (100 * tol))
    # The physics auditor runs on the photomask file (opencif.checks.auditor). Its hard checks are a pass condition in
    # every mode, and its M14 (the Reynolds and Stokes gates on the drawn mask) is part of the validity envelope.
    pa = V.get('physics_audit')
    if pa is None:
        chk('physics audit of the photomask file (P0-P10)', False, 'not run', 'the audit must run', 'incomplete')
    else:
        drc_missing = not V.get('drc') or 'error' in V.get('drc', {})
        chk('physics audit of the photomask file (P0-P10)', not pa['hard_failures'],
            dict(status=pa['status'], failed=pa['hard_failures'],
                 checks={c['id']: ('ok' if c['ok'] else ('reported' if c['role'] == 'report' else 'FAIL'))
                         for c in pa['checks']}),
            'the photomask read back from its file passes the physics auditor\'s hard checks P0, R0, P1-P6b, P9, R6b '
            'and P10, and P9-AR in strict mode (CIF_v2\'s P1-P10 ported; P0: the file is the drawing; R0: the '
            'auditor\'s own ruler on the straight rows; R6b: the posts; P7 and P8 are reported, M14 is part of the '
            'validity envelope; validation.physics_audit has each check, its value and its requirement)',
            'incomplete' if (drc_missing and pa['hard_failures'] == ['P9']) else 'hard')
        p8_ = next(c for c in pa['checks'] if c['id'] == 'M14')
        chk('physics audit M14: the design code\'s Reynolds and Stokes gates on the drawn mask', p8_['ok'],
            dict((k, p8_['value'][k]) for k in ('Re_channel_max', 'Re_gap_max', 'Re_delta_max')), p8_['requirement'],
            'envelope')
    dcr = V.get('design_code_reports')
    if dcr:                                               # reported, as in the design code
        vo_ = dcr['V_operating_regime']['at_5_uL_min']
        chk('the design code\'s reports and acceptance validators (reported, not a pass condition)', True,
            dict(dP_device_kPa=dcr['B3_pressure_drop']['drawn_device_kPa'], eta_F=dcr['B6_F_real']['eta_F'],
                 transfer_bias=dcr['B7_transfer_bias']['D_true_over_D_predicted'],
                 V_operating_regime_at_5_uL_min=vo_['ok'], V_steric=dcr['V_steric']['ok'],
                 V_fab_tolerance_yield=dcr['V_fab_tolerance']['yield_F']),
            'design()\'s mandatory reports (pressure drop, eta F, transfer bias, sensitivity map, wall porosity) and '
            'run_acceptance\'s operating validators (V_operating_regime at 5 uL/min, V_steric, V_fab_tolerance), each by '
            'the design code\'s formula: validation.design_code_reports', 'model')
    if rc['declared_rows']:
        chk('declared departures (reported, not a pass condition)', True,
            dict(by_option=rc.get('declared_by_option'), reverse_flow_share_of_inlet=round(rc['reverse_flow_share_of_inlet'], 5)),
            'rows each declared option moves, inside its zone, predicted by the design network before measuring', 'declared')
    if chip_rect:
        # The packaging budget is reported and is never a pass condition. App. B's thesis says that "it reports, without
        # refusing, a design over its length budget (C8, CQ5)", and CQ5 says that the budget check "reports the result
        # without failing".
        b = fluid.bounds
        cw, ch = chip_rect
        fx, fy = b[2] - b[0], b[3] - b[1]
        chk('fits the chip area (reported, not a pass condition: App. B thesis; C8 CQ5)',
            fx <= cw - 2 * opts['margin'] + 1e-6 and fy <= ch - 2 * opts['margin'] + 1e-6,
            [round(fx / 1e3, 3), round(fy / 1e3, 3)], 'inside %.2f x %.2f mm less a %.0f um margin' % (cw / 1e3, ch / 1e3, opts['margin']),
            'budget')
    # The aspect ratio of the narrowest drawn channel is checked against the process limit (CIF_v2
    # w_min = max(F_min, h/AR_max); spec v5 G7 "widths >= w_min", and G11, which folds a side channel narrower than
    # w_min, a fold that opencif does not draw). In strict mode G7 is a pass condition on the drawing, and the planned
    # straight side channels are already refused before drawing by wmin_check. Evaluate and reproduce modes report it
    # (docs/LIMITATIONS.md 7). The inventory is the drawn channels' measured widths, which are the gap openings
    # (measured, linear rows), the center channel (measured at both ends of its taper), every straight side section, the
    # U-turn sides, the head units, and the outlet leads. The inlet stem and filter are checked by the DRC at the floor
    # only (reviews of rc4).
    widths = AU.drawn_widths(V, spec)
    g_ = V.get('gates') or {}
    ua_ = g_.get('plateau_u_min_all_planned_rows')
    pre_ = 'declared rows below the kappa plateau'
    if ua_ is not None and g_.get('plateau_ok') and ua_ < S.U_LO:      # declared rows below the kappa plateau
        _renote(notes, pre_, '%s: planned row cutoffs down to u = D_row/w_c = %.4f < %.3f (the cutoff law is not '
                             'calibrated there; %.1f %% of the extracted flow at such rows); the device\'s own cutoff D is '
                             'on the plateau' % (pre_, ua_, S.U_LO, 100.0 * float(g_.get('plateau_F_share_below_all_planned_rows') or 0.0)))
    else:
        _renote(notes, pre_, None)
    if opts.get('passage_report', True):
        rcm_ = rc.get('D_realised_max_um') if isinstance(rc, dict) else None
        # The passage comparison is a report rather than a check (ADJ-C9, R3).
        V['passage'] = passage_report(widths, S.KAPPA_G * spec.delta_star, notes, D_max=rcm_,
                                      G=widths.get('gap opening G'),
                                      rows_no_cutoff=(rc.get('rows_f_half_or_more', 0) if isinstance(rc, dict) else 0))
    w_where = min(widths, key=widths.get)
    ar_drawn = spec.h / widths[w_where]
    ar_max = opts.get('AR_max', BU.AR_MAX_DEFAULT)
    strict = opts['mode'] == 'strict'
    chk('aspect ratio of the narrowest CIF channel (gap openings, centre channel, side channels, U-turn sides, head '
        'units, outlet leads, measured; %s: '
        'LIMITATIONS 7)' % ('a pass condition in strict mode' if strict else 'reported, not a pass condition in %s mode'
                            % opts['mode']),
        ar_drawn <= ar_max + 1e-9, dict(depth_over_width=round(ar_drawn, 2), at=w_where, width_um=round(widths[w_where], 3)),
        'depth/width <= process.aspect_ratio_max = %g (CIF_v2 w_min = max(F_min, h/AR_max) = %.2f um; spec v5 G7, G11)'
        % (ar_max, max(spec.F_min, spec.h / ar_max)), 'hard' if strict else 'fabrication')
    g = V['gates']
    env = dict(plateau=g['plateau_ok'], M13=g['M13_ok'], Reynolds=g['Re_ok'], Stokes=g['stokes_ok'],
               foot=g['plateau_u_min'] >= S.U_FF, lamina=g.get('lamina_ok', True) and g.get('lamina_all_rows_ok', True),
               taper=g.get('taper_ok', True))
    chk('validity envelope (kappa plateau, gap dominance, Reynolds, Stokes gap and lamina, thin lamina, taper rate)',
        all(env.values()),
        dict(u_min=round(g['plateau_u_min'], 4), M13=round(g['M13_gap_dominance_max'], 4), Re_max=round(g['Re_max'], 3),
             Re_gap_max=float('%.3g' % g['Re_gap_max']), Re_delta_max=float('%.3g' % g['Re_delta_max']),
             lamina_over_depth=float('%.3g' % g['lamina_over_depth']),
             lamina_over_depth_planned_max=float('%.3g' % g['lamina_over_depth_planned_max']),
             taper_rate_max=float('%.3g' % g.get('taper_rate_max', 0.0)),
             f_design_max=float('%.3g' % g['f_design_max']), cutoff_band_um=g['cutoff_band_um']),
        'u_min >= %.3f (no gap below the kappa knee); M13 < %g; Re <= %g; Re_gap and Re_delta <= %g (CIF_v2 M13/M14); '
        'delta(i)/h <= %g at every planned row (E11); |w_c(i) - w_c(i-1)|/(G + g_p) <= %g (M16)'
        % (S.U_LO, BU.M13_MAX, BU.RE_MAX, BU.RE_STOKES, BU.EPS_LAM, BU.EPS_TAPER), 'envelope')
    if (g.get('Re_gap_planned_max', 0) > BU.RE_STOKES or g.get('Re_delta_planned_max', 0) > BU.RE_STOKES) and g['stokes_ok'] \
            and not any('planned away from f*' in n for n in notes):
        notes.append('rows planned away from f* (declared departures, e.g. an unsealed gap 1) reach Re_gap %.2e and '
                     'Re_delta %.2e, above the Stokes gate %g: reported, not gated, since those rows claim no cutoff at D '
                     '(CIF_v2 M14 gates cutoff rows only)' % (g['Re_gap_planned_max'], g['Re_delta_planned_max'], BU.RE_STOKES))
    if g.get('f_planned_max', 0) > S.FGAP_CAP and not g.get('f_beyond_calibration') \
            and not any('planned f_gap' in n for n in notes):
        notes.append('a row planned away from f* (a declared departure) has planned f_gap %.2e, above kappa_g\'s calibration '
                     'cap %.2e: its reported cutoff is extrapolated' % (g['f_planned_max'], S.FGAP_CAP))
    if g.get('f_beyond_calibration') and not any('calibration cap' in n for n in notes):
        notes.append('the largest design f_gap, %.2e, is above kappa_g\'s calibration cap %.2e (glossary \\fgapDesignCap, '
                     'Gifford\'s array): the cutoff law is extrapolated on those rows (reported, not gated)'
                     % (g['f_design_max'], g['f_calibration_cap']))
    if g.get('f_above_soft') and not any('f_soft' in n for n in notes):
        notes.append('the largest design f_gap, %.2e, is above spec v5\'s soft bound f_soft = %.0e (CIF_v2 Pairing; spec v5 '
                     'sec. 18 and sec. 19 D3 "WARN": outside the best-justified regime; reported, not gated)' % (g['f_design_max'], g['f_soft']))
    ms = V.get('model_scenarios') or {}
    model_txt = []
    if ms and 'error' not in ms:
        val = {}
        b1, b2 = ms.get('junction_band_lo'), ms.get('junction_band_hi')
        if b1 and b2:
            nb = max(b1['rows_newly_outside_tol'], b2['rows_newly_outside_tol'])
            nl = max(b1.get('rows_losing_their_cutoff', 0), b2.get('rows_losing_their_cutoff', 0))
            mx = max(b1['max_cutoff_change'] or 0.0, b2['max_cutoff_change'] or 0.0)
            rows_ = sorted(set(b1['newly_outside_rows']) | set(b2['newly_outside_rows']))
            val['junction band (bend, head T-junction and outlet stretch-end terms at both ends of their range)'] = dict(
                rows_newly_outside_tol=[b1['rows_newly_outside_tol'], b2['rows_newly_outside_tol']],
                rows_losing_their_cutoff=nl, max_cutoff_change_over_D=round(mx, 4), rows=_rows_txt(rows_, None))
            for lab in ('t_junction_ignored', 'junctions_ignored', 'bend_correction_off'):
                if ms.get(lab):
                    val[lab] = dict(rows_newly_outside_tol=ms[lab]['rows_newly_outside_tol'],
                                    rows_losing_their_cutoff=ms[lab].get('rows_losing_their_cutoff', 0),
                                    max_cutoff_change_over_D=round(ms[lab]['max_cutoff_change'] or 0.0, 4))
            sh = max(b1.get('share_of_F_in_those_rows') or 0.0, b2.get('share_of_F_in_those_rows') or 0.0)
            val['junction band (bend, head T-junction and outlet stretch-end terms at both ends of their range)']['share_of_F_in_those_rows'] = sh
            if nb:
                model_txt.append('junction band: %d row(s) beyond +-%g%%%s (up to %.1f%% of D; rows %s; %.2g%% of F)'
                                 % (nb, 100 * tol, (', %d losing their cutoff' % nl) if nl else '', 100 * mx,
                                    _rows_txt(rows_), 100 * sh))
        elif ms.get('bend_model', {}).get('model') == 'none':
            model_txt.append('plain 1-D bends (bend_model: none): comparison only')
        ob = [ms.get('outlet_side_resistance_+1%'), ms.get('outlet_side_resistance_-1%')]
        if all(ob):
            val['outlet balance (both side outlets +-1 %)'] = dict(
                rows_newly_outside_tol=[ob[0]['rows_newly_outside_tol'], ob[1]['rows_newly_outside_tol']],
                rows_losing_their_cutoff=max(ob[0].get('rows_losing_their_cutoff', 0), ob[1].get('rows_losing_their_cutoff', 0)),
                max_cutoff_change_over_D=round(max(ob[0]['max_cutoff_change'] or 0.0, ob[1]['max_cutoff_change'] or 0.0), 4),
                rows=_rows_txt(sorted(set(ob[0]['newly_outside_rows']) | set(ob[1]['newly_outside_rows'])), None),
                share_of_F_in_those_rows=max(ob[0].get('share_of_F_in_those_rows') or 0.0,
                                             ob[1].get('share_of_F_in_those_rows') or 0.0),
                dDlast_dlnRside_out=round(ms.get('dDlast_dlnRside_out', 0.0), 4))
        tub = {k: v for k, v in ms.items() if k.startswith('tubing_')}
        if tub:
            nt = max(v['rows_newly_outside_tol'] for v in tub.values())
            val['outlet tubing'] = dict({k: [v['rows_newly_outside_tol'], round(v['max_cutoff_change'] or 0.0, 4)]
                                         for k, v in tub.items()}, dDlast_dlnRextc=round(ms.get('dDlast_dlnRextc', 0.0), 4))
            if nt:
                model_txt.append('last rows depend on the outlet tubing (%d rows beyond +-%g%% in its scenarios)' % (nt, 100 * tol))
        chk('1-D model band (reported, not a pass condition)', True, val,
            'cutoff change on the drawn network when the junction terms (bends, head T-junctions, the outlet stretch\'s '
            'end) take the ends of their range (the depth weight moved; the T-junction term widened; the stretch-end term '
            '+-30 %) (and, for information, ignored, and with no bend correction), when both side outlets\' resistance '
            'moves +-1 % (the outlet balance the last rows rest on), and when the off-chip tubing varies (length '
            'x0.5/x2, centre I.D. +-0.025 mm); counted on rows inside D +- tol on the drawn mask', 'model')
    elif ms:
        chk('1-D model band (reported, not a pass condition)', False, ms.get('error'), 'the scenarios must run', 'model')
    hard_fail = [c['check'] for c in out if not c['ok'] and c['kind'] == 'hard']
    incomplete = [c['check'] for c in out if not c['ok'] and c['kind'] == 'incomplete']
    env_fail = [c['check'] for c in out if not c['ok'] and c['kind'] == 'envelope']
    over_budget = [c for c in out if not c['ok'] and c['kind'] == 'budget']
    over_ar = [c for c in out if not c['ok'] and c['check'].startswith('aspect ratio')]
    band = g['cutoff_band_um']
    why = []
    if not env['plateau'] or not env['foot']:
        why.append('u_min %.4f below the kappa plateau %.3f' % (g['plateau_u_min'], S.U_LO))
    if not env['M13']:
        why.append('gap dominance M13 %.3f >= 0.05' % g['M13_gap_dominance_max'])
    if not env['Reynolds']:
        why.append('Re %.1f > 10: inertial (Dean) flow at bends and gap mouths is outside the Stokes model' % g['Re_max'])
    if not env['Stokes']:
        why.append('Stokes gates: Re_gap %.2e, Re_delta %.2e > %g: the skimming junction is not Stokes'
                   % (g['Re_gap_max'], g['Re_delta_max'], BU.RE_STOKES))
    if not env['lamina']:
        why.append('thin lamina E11: delta/h %.4f > %g: the near-wall quadratic the cutoff law rests on is not valid'
                   % (max(g['lamina_over_depth'], g['lamina_over_depth_planned_max']), BU.EPS_LAM))
    if not env['taper']:
        why.append('taper rate M16: |dw_c|/(G + g_p) %.3e at row %d > %g: the flow is not quasi-developed'
                   % (g['taper_rate_max'], g['taper_rate_row'], BU.EPS_TAPER))
    if pa is not None and env['Reynolds'] and env['Stokes'] and not p8_['ok']:
        why.append('on the drawn mask (physics audit M14): Re %.3g, Re_gap %.2e, Re_delta %.2e (limits %g, %g)'
                   % (p8_['value']['Re_channel_max'], p8_['value']['Re_gap_max'], p8_['value']['Re_delta_max'],
                      BU.RE_MAX, BU.RE_STOKES))
    env_txt = ('OUTSIDE THE VALIDITY ENVELOPE (%s): the cutoff law is not calibrated there; %s'
               % ('; '.join(why), band_txt(band)))
    if hard_fail:
        verdict, code = 'FAIL (%s)' % '; '.join(hard_fail), 2
    elif incomplete:
        verdict, code = 'INCOMPLETE (%s)' % '; '.join(incomplete), 4
    else:
        dec_names = [k for k, v in rc['declared_rows'].items() if v]
        if opts['mode'] == 'reproduce':
            head = 'REPRODUCED' + (' with declared departures (%s)' % ', '.join(dec_names) if dec_names else '')
        else:
            head = 'PASS' + (' with an allowed departure (%s)' % ', '.join(dec_names) if dec_names else '')
        verdict, code = (head + ' -- ' + env_txt, 3) if env_fail else (head, 0)
    # An outlet stretch-end beyond the solved 2-D tables is refused in strict mode, and here it is held and reported.
    td_ = (((rep.get('design') or {}).get('outlet') or {}).get('stretch_end') or {}).get('table_domain')
    if td_ and td_ != 'inside the solved tables':
        model_txt.append('stretch-end term outside the solved tables: %s' % td_)
    if model_txt and code in (0, 3):
        verdict += ' [1-D model: %s]' % '; '.join(model_txt)
    if over_budget:
        verdict += ' -- OVER THE CHIP AREA (%s mm drawn; reported, not refused or failed)' % '; '.join(
            '%.2f x %.2f' % tuple(c['value']) for c in over_budget)
    if over_ar and not strict:
        verdict += (' -- ABOVE THE PROCESS ASPECT-RATIO LIMIT (depth/width %.1f > %g at the narrowest %s; reported, not '
                    'refused or failed in %s mode)' % (over_ar[0]['value']['depth_over_width'], ar_max,
                                                       over_ar[0]['value']['at'], opts['mode']))
    status_ar = ('ABOVE THE LIMIT' + (' (a FAIL in strict mode)' if strict else ' (reported)')) if over_ar else 'WITHIN THE LIMIT'
    # Spec v5 G9 compares each leg's path length with the chip's leg length, and the result is reported without
    # failing (App. B).
    lb = V.get('length_budget') or {}
    if chip_rect and lb.get('leg_path_lengths_mm'):
        ell = (chip_rect[0] - 2 * opts['margin']) / 1e3
        lb['leg_budget_mm'] = round(ell, 3)
        lb['legs_over_budget'] = [j + 1 for j, x in enumerate(lb['leg_path_lengths_mm']) if x > ell + 1e-9]
        status_g9 = 'OVER (legs %s; reported)' % lb['legs_over_budget'] if lb['legs_over_budget'] else 'WITHIN'
    else:
        status_g9 = 'NOT GIVEN (no chip)'
    # The verdict's parts are kept as separate fields, so that no single word stands for all of them.
    status = dict(hard_checks=('FAIL' if hard_fail else ('INCOMPLETE' if incomplete else 'PASS')),
                  validity_envelope=('OUTSIDE' if env_fail else 'INSIDE'),
                  chip_area=('NOT GIVEN' if not chip_rect else ('OVER' if over_budget else 'FITS')),
                  aspect_ratio_of_CIF_channels=status_ar,
                  length_budget_G9=status_g9,
                  drawn_mask_assessment='RUN: opencif\'s own checks on the drawn polygons (DRC at the floor, topology, '
                                        'continuity, ports, clearances, corner radii, aspect ratio) and the 1-D network '
                                        'measured on them (Q1, Q2)',
                  physics_audit=_audit_status(V.get('physics_audit'), V.get('drc')),
                  fabricated_or_measured='NO: Q1 and Q2 are statements of the 1-D network on the drawn mask')
    if opts['mode'] == 'reproduce':
        tf_ = opts.get('targets_fitted')          # the targets that the spec fitted to the paper's figures
        status['reproduction'] = ('DESIGN-LEVEL: %s, drawn with the paper\'s parts and printed dimensions where printed; '
                                  'not a full-mask match, not a measured-performance match (docs/REPRODUCTIONS.md)' % (
                                      ('targets fitted to the paper\'s figures where they contradict its printed values '
                                       '(paper.targets_fitted: %s), the others from printed values' % tf_) if tf_ else
                                      'D and F from printed values'))
    else:                                                 # one status field per status code (R16)
        status['reproduction'] = 'NOT CLAIMED: not a reproduction (gates.mode is not reproduce)'
    # The status carries a pointer to the trace of this device against its source figures (the spec's
    # paper.source_trace), which the pipeline neither reads nor checks, and every status field as one enumerated code
    # (STATUS_CODES) (ADJ-C6).
    st_ = opts.get('source_trace')
    # The code is NAMED rather than TRACED, because the build cannot know the trace, which is made from its own
    # photomask afterwards. tests/test_rc7.py checks that each shipped preset's pointer resolves to traces of its
    # shipped mask (codereview2 N1).
    status['source_trace'] = (('NAMED: %s (a pointer given by the spec, paper.source_trace; the build neither reads '
                               'nor checks it)' % st_) if st_ else
                              'NOT NAMED: the spec names no trace of this device against its source figures')
    codes = dict(hard_checks=status['hard_checks'], validity_envelope=status['validity_envelope'],
                 chip_area=status['chip_area'].replace(' ', '_'),
                 aspect_ratio_of_CIF_channels=('ABOVE_LIMIT' if over_ar else 'WITHIN_LIMIT'),
                 length_budget_G9=('NOT_GIVEN' if status_g9.startswith('NOT GIVEN') else
                                   ('OVER' if status_g9.startswith('OVER') else 'WITHIN')),
                 drawn_mask_assessment='RUN', physics_audit=_audit_code(V.get('physics_audit'), V.get('drc')),
                 fabricated_or_measured='NO',
                 reproduction=('DESIGN_LEVEL' if opts['mode'] == 'reproduce' else 'NOT_CLAIMED'),
                 source_trace=('NAMED' if st_ else 'NOT_NAMED'))
    for k_, v_ in codes.items():                          # an explicit check, because -O strips an assert (R16)
        if v_ not in STATUS_CODES[k_]:
            raise ValueError('status code %s = %r is not one of %s' % (k_, v_, STATUS_CODES[k_]))
    return dict(verdict=verdict, exit_code=code, checks=out, notes=notes, cutoff_tolerance=tol, mode=opts['mode'],
                status=status, status_codes=codes, over_chip_area=bool(over_budget), over_aspect_ratio=bool(over_ar),
                reproduction_only_options=opts['repro_only'], defaults_used=opts.get('defaults_used', []))


def summary_md(rec, path):
    """Write the record's Markdown summary to path, with the verdict, the targets against the drawn mask, the design,
    every check, the physics audit, the status, and the notes."""
    a = rec['assessment']
    d = rec['design']
    V = rec['validation']
    rc = V['realised_cutoff']
    lines = ['# %s' % d['name'], '', '**Verdict: %s**' % a['verdict'], '',
             '| | target | on the drawn mask |', '|---|---|---|',
             '| cutoff D (µm) | %.4f | %.4f – %.4f over all gaps |' % (rc['D_target_um'], rc['D_realised_min_um'],
                                                               rc['D_realised_max_um']),
             '| gaps outside D ± %g%% (not declared) | 0 | %d |' % (100 * a['cutoff_tolerance'], rc['undeclared_rows_outside']),
             '| share of the net filtrate extracted within D ± %g%% | | %.3f%% |' % (
                 100 * a['cutoff_tolerance'], 100 * rc['share_of_F_within_%g%%' % (100 * a['cutoff_tolerance'])]),
             '| filtrate fraction F | %.5f | %.6f |' % (V['network']['F_target'], V['network']['F_as_drawn']),
             '', 'Design: N = %d gaps (%d in the meander head), layout %s with %d leg(s), footprint %.2f × %.2f mm.' % (
                 d['N'], d['NM'], d['layout'], d['n_legs'], V['drawn_parameters']['footprint_mm'][0], V['drawn_parameters']['footprint_mm'][1]),
             '', 'Cutoff law: D = κ_g δ*, κ_g = %.4f, 95 %% CI [%.3f, %.3f] (Gifford 2014 calibration; App. B C7); %s.' % (
                 BU.KERNEL_PAIR['kappa_g'], BU.KERNEL_PAIR['kappa_ci'][0], BU.KERNEL_PAIR['kappa_ci'][1],
                 band_txt(V['gates']['cutoff_band_um'])),
             '', '| check | result | value | requirement |', '|---|---|---|---|']
    for c in a['checks']:
        lines.append('| %s | %s | `%s` | %s |' % (c['check'], ('ok' if c['ok'] and c['kind'] not in ('declared', 'model') else
                                                              ('reported' if c['kind'] in ('declared', 'model') else
                                                               ('**over (reported)**' if c['kind'] in ('budget', 'fabrication') else
                                                                ('**fail**' if c['kind'] == 'hard' else c['kind'])))),
                                               json.dumps(c['value'], default=float), c['requirement']))
    pa = V.get('physics_audit')
    if pa:                                                 # the physics audit, check by check
        lines += ['', '**Physics audit of %s: %s** (opencif.checks.auditor; each check\'s value and requirement are in the '
                  'record, validation.physics_audit)' % (pa['photomask'].get('file_final', pa['photomask']['file']), pa['status']), '',
                  '| check | | result | role |', '|---|---|---|---|']
        for c in pa['checks']:
            lines.append('| %s | %s | %s | %s |' % (c['id'], c['check'], 'ok' if c['ok'] else
                                                    ('reported' if c['role'] == 'report' else '**fail**'), c['role']))
    if a.get('status_codes'):                              # the enumerated status codes (ADJ-C6)
        lines += ['', 'Status codes: ' + ', '.join('%s %s' % (k, v) for k, v in a['status_codes'].items()) + '.']
    if a.get('status'):
        lines += ['', 'Status: ' + '; '.join('%s %s' % (k.replace('_', ' '), v) for k, v in a['status'].items()) + '.']
    if a['notes']:
        lines += ['', 'Notes:'] + ['- ' + n for n in a['notes']]
    if a.get('defaults_used'):
        lines += ['', 'Values not given in the spec (defaults and derivations):'] + ['- ' + n for n in a['defaults_used']]
    open(path, 'w').write('\n'.join(lines) + '\n')


def prepare(d, say=lambda m: None):
    """Do everything that comes before drawing, which `build` and `check` share, and return
    (spec, opts, sc, gates, envelope_failures, notes). The steps are the mapping of the spec onto the engine, the seal
    policy, the head width, the floor start, the floor guard without a head, the design plan, and its gates."""
    spec, opts = SP.to_engine(d)
    L.WALL_RULE = opts['wall_rule']
    try:
        return _prepare(spec, opts, say)
    except SP.SpecError:
        raise
    # A plan that cannot be completed is a spec error rather than a crash, whether the closure has no solution
    # (ValueError) or the march needs more than its 400,000 rows (RuntimeError).
    except (ValueError, RuntimeError) as e:
        fs = S.f_of(spec.delta_star, spec.wc0)
        raise SP.SpecError('the design cannot be planned (D %.3g um on a %.3g um centre channel: u = D/w_c = %.3f, '
                           'f_gap %.3g at gap 1): %s: %s%s%s' % (
                               S.KAPPA_G * spec.delta_star, spec.wc0, S.KAPPA_G * spec.delta_star / spec.wc0, fs,
                               type(e).__name__, e,
                               ('.  f_gap is %.0f x kappa_g\'s calibration cap %.2e: reduce D or widen the centre channel'
                                % (fs / S.FGAP_CAP, S.FGAP_CAP)) if fs > S.FGAP_CAP else '',
                               ('.  Also found before planning: %s' % opts['preflight_refusal'])
                               if opts.get('preflight_refusal') else ''))


def gamma_gs_check(spec, opts, notes):
    """Check the fold width GS (meander channel_um) over the gap G against gamma_GS, which is the pre-flight
    "ASSERT GS/G in gamma_GS" of spec v5 sec. 18 (and sec. 14.4; in CIF_v2, design() raises GateError).

    The ratio is checked only when a meander head is drawn, because without a head GS is not drawn. Strict mode refuses
    a ratio outside the band (opts['preflight_refusal']), and evaluate and reproduce modes report it in the notes. A
    narrow head's narrow units are not GS-wide (spec v5 sec. 14.2 draws arms at GS), which the notes state without
    gating it."""
    if spec.head == 'none':
        return
    lo, hi = BU.GAMMA_GS
    ratio = spec.GS / spec.G
    if not (lo - 1e-9 <= ratio <= hi + 1e-9):
        msg = ('fold width over gap GS/G = %.1f/%.1f = %.3f, outside gamma_GS [%g, %g] (meander.channel_um over '
               'process.gap_um; spec v5 sec. 14.4 and its sec. 18 pre-flight gate, CIF_v2 design(); glossary: a supplied process '
               'parameter)' % (spec.GS, spec.G, ratio, lo, hi))
        _preflight(opts, notes, msg)
    if getattr(spec, 'w_narrow', None):
        notes.append('narrow head: the units that would need GS-wide arms are drawn %.1f um wide (w/G = %.2f); spec v5 '
                     'sec. 14.2 draws arms at GS, and gamma_GS [%g, %g] is checked on GS only (GS/G = %.3f)'
                     % (spec.w_narrow, spec.w_narrow / spec.G, lo, hi, ratio))


def _preflight(opts, notes, msg):
    """Apply a pre-flight gate of spec v5 sec. 18 (CIF_v2 design()). Strict mode refuses the design
    (opts['preflight_refusal'], with several refusals joined), and evaluate and reproduce modes report the gate in the
    notes."""
    if opts['mode'] == 'strict':
        opts['preflight_refusal'] = '; '.join([m for m in (opts.get('preflight_refusal'), msg) if m])
    else:
        notes.append(msg + ': reported, not refused (gates.mode = %s)' % opts['mode'])


def wmin_check(spec, opts, notes, sc=None):
    """Check CIF_v2's aspect-ratio floor w_min = max(F_min, h/AR_max) (Process.w_min, "S12: the smallest buildable
    channel width"; spec v5 sec. 18) where CIF_v2 enforces it.

    Before planning (sc None), the check is the pre-flight "GS or G below w_min", with GS checked only when a meander
    head is drawn, because without one GS is not drawn. Once the design is planned, the checks are the terminal side
    width w_s(N) >= w_min (design() raises GateError) and every center width w_c >= w_min (its validator G7), together
    with the planned straight side channels. Strict mode refuses, and the other modes report. AR_max is
    process.aspect_ratio_max (by default CIF_v2's 10). On the drawn mask, assess() checks the narrowest measured
    channel's aspect ratio, which is a pass condition in strict mode and is reported in the other modes
    (docs/LIMITATIONS.md 7)."""
    w_min, ar = opts['w_min'], opts['AR_max']
    src = ('w_min = max(F_min, h/AR_max) = max(%.2f, %.1f/%g) = %.2f um (CIF_v2 Process.w_min and design(); spec v5 sec. 18; '
           'process.aspect_ratio_max)' % (spec.F_min, spec.h, ar, w_min))
    if sc is None:
        if spec.G < w_min - 1e-9:
            _preflight(opts, notes, 'gap G = %.2f um below %s' % (spec.G, src))
        if spec.head != 'none' and spec.GS < w_min - 1e-9:
            _preflight(opts, notes, 'meander channel GS = %.2f um below %s' % (spec.GS, src))
        return
    ws = np.asarray(sc.ws, float)
    if np.any(np.isfinite(ws)):
        wsN = float(ws[np.isfinite(ws)][-1])
        if wsN < w_min - 1e-9:
            _preflight(opts, notes, 'terminal side width w_s(N) = %.2f um below %s' % (wsN, src))
    wc_min = float(np.min(np.asarray(sc.wc, float)))                     # CIF_v2's validator G7, every w_c >= w_min
    if wc_min < w_min - 1e-9:
        _preflight(opts, notes, 'centre width w_c = %.2f um below %s' % (wc_min, src))
    # For the straight side channels, spec v5 G7 asks "widths >= w_min" of the drawing, and G11 draws a side channel
    # narrower than w_min as a fold, which opencif does not draw. Therefore strict mode refuses a design whose planned
    # straight side channel falls below w_min, and the other modes report it. The drawn mask's aspect-ratio check in
    # assess() is then a pass condition in strict mode and a report in the others.
    nm, n_ = int(getattr(sc, 'NM', 0)), int(getattr(sc, 'N', len(ws)))
    ws_str = ws[nm:n_ - 1]
    ws_str = ws_str[np.isfinite(ws_str)]
    if ws_str.size and float(np.min(ws_str)) < w_min - 1e-9:
        k = int(nm + np.argmin(np.where(np.isfinite(ws[nm:n_ - 1]), ws[nm:n_ - 1], np.inf)))
        _preflight(opts, notes, 'straight side channel w_s = %.2f um (row %d) below %s: spec v5 G7 asks widths >= w_min, and '
                   'G11 would draw that row as a fold, which opencif does not draw' % (float(np.min(ws_str)), k + 1, src))


def refusal(opts, env_bad):
    """Apply the one refusal rule that run(), `opencif check`, and generate(check_only=True) share, and return
    (refused, reasons, tag). A seal refusal, an outlet refusal, and a layout refusal (parts that collide on the outline)
    refuse in any mode, whereas the pre-flight gates and the validity envelope refuse in strict mode."""
    strict_env = bool(env_bad) and opts['mode'] == 'strict'
    reasons = (([opts['seal_refusal']] if opts.get('seal_refusal') else []) +
               ([opts['outlet_refusal']] if opts.get('outlet_refusal') else []) +
               ([opts['layout_refusal']] if opts.get('layout_refusal') else []) +
               ([opts['preflight_refusal']] if opts.get('preflight_refusal') else []) + (list(env_bad) if strict_env else []))
    tag = ' and '.join((['seal'] if opts.get('seal_refusal') else []) + (['outlet'] if opts.get('outlet_refusal') else []) +
                       (['layout'] if opts.get('layout_refusal') else []) +
                       (['strict'] if (strict_env or opts.get('preflight_refusal')) else []))
    return bool(reasons), reasons, tag


def layout_remedy(spec, sc, kind, n_legs, opts, coll):
    """Say what would separate the colliding parts, for the refusal's advice (the clearance guarantee).

    The advice names the parts and the shortfall, the layout values that the spec gave (which the generator then uses
    as given), and, for a paperclip, the leg spacing that the generator's own rule gives this design
    (layout._plan_paperclip), which is

        max(2 hw_max, hw_head + hw_max) + 150 um

    where hw_max is the widest linear row's half-width and hw_head is the head's."""
    out = []
    pc = coll['part_clearance']
    if not pc['ok']:
        c0 = pc['closest'][0]
        out.append("move '%s' and '%s' apart until at least the %.1f um floor separates them (%s)"
                   % (c0['parts'][0], c0['parts'][1], pc['floor_um'],
                      'they overlap now' if c0['gap_um'] <= 0.0 else 'now %.1f um' % c0['gap_um']))
    given = []
    if opts.get('leg_spacing') is not None:
        given.append('layout.leg_spacing_um = %g' % float(opts['leg_spacing']))
    if getattr(spec, 'leg_Y', None) is not None:
        given.append('reproduction.leg_Y = %s' % [round(float(v), 1) for v in spec.leg_Y])
    if getattr(spec, 'pc_inlet_offset', 0.0):
        given.append('reproduction.pc_inlet_offset = %g' % float(spec.pc_inlet_offset))
    if given:
        out.append('the layout uses the spec\'s %s as given' % ', '.join(given))
    pair_ = set(pc['closest'][0]['parts']) if not pc['ok'] else set()
    ext_ = (coll.get('extensions') or {}).get('tried') or []
    if ext_:
        out.append('the outlet stretch was lengthened up to %.0f um and the outline still collides' % max(
            t.get('extension_um', 0.0) for t in ext_))
    if (not coll['port_spacing']['ok'] or not coll['separation']['separated'] or 'outlet' in pair_) and \
            kind == 'paperclip' and n_legs > 1 and n_legs % 2 == 0:
        out.append('the last leg ends beside the inlet: an odd leg count ends a paperclip at the far end (legs: auto '
                   'skips the counts whose parts collide)')
    if (not coll['port_spacing']['ok'] or 'outlet' in pair_) and (kind == 'linear' or n_legs <= 1):
        out.append('a longer outlet stretch (outlet.stretch_um, now %.0f um) moves the outlet ports away from the inlet '
                   'port' % float(getattr(spec, 'L_out', 0.0) or 0.0))
    if given and kind == 'paperclip' and n_legs > 1 and sc is not None:
        N, NM = sc.N, sc.NM
        hw_max = float(np.nanmax(sc.wc[NM:N - 1] / 2.0 + spec.l_g + sc.ws[NM:N - 1]))
        hw_head = float(L.head_halfwidth(sc))
        d_auto = max(2.0 * hw_max, hw_head + hw_max) + 150.0
        out.append('the generator\'s own rule spaces this paperclip\'s legs %.0f um apart (max(2 x %.0f, %.0f + %.0f) + '
                   '150 um: the widest row\'s half-width and the head\'s); give no spacing to use it, or space the '
                   'legs beside the head at least that far' % (d_auto, hw_max, hw_head, hw_max))
    out.append(('or ' if out else '') + 'change the parts themselves (a narrower head, smaller ports, another leg '
               'count or layout)')
    return '; '.join(out)


def layout_check(spec, sc, opts, notes, say=lambda m: None):
    """Make in memory the layout and outline that a build would draw, with nothing written, and apply the collision
    checks that a build applies before it writes the photomask (build.collisions; the clearance guarantee).

    `opencif check` and generate(check_only=True) use it, so that a check refuses exactly what a build would refuse. It
    sets opts['layout_refusal'] (or opts['outlet_refusal'] when the outlet turns out not to close on the fitted layout)
    and returns the layout summary for the check's record."""
    import copy
    sp = copy.deepcopy(spec)
    try:
        kind, n_legs, sp, sc_ = fit_layout(sp, sc, opts, notes, say)
        sp = copy.deepcopy(sp)
        fluid, sc_f, info, _, _, _, coll = BU.plan_and_draw(sp, kind, n_legs, opts['leg_spacing'], opts['failsafe'],
                                                           say, 'check', time.time())
    except L.OutletOutsideTables as e:
        opts['outlet_refusal'] = ('the %s outlet: %s: build it in evaluate or reproduce mode, where the end values '
                                  'are held and reported, or change the outlet' % (opts.get('outlet'), e))
        return None
    except L.OutletInfeasible as e:
        opts['outlet_refusal'] = 'the %s outlet cannot close this design (F %.4f): %s' % (opts.get('outlet'),
                                                                                          spec.F_target, e)
        return None
    b = fluid.bounds
    out = dict(layout=kind, legs=n_legs, footprint_mm=[(b[2] - b[0]) / 1e3, (b[3] - b[1]) / 1e3],
               leg_Y=sc_f.info.get('leg_Y'), leg_d=sc_f.info.get('leg_d'), R_turn=sc_f.info.get('R_turn'),
               outlet_stretch_extended_um=float(getattr(sp, 'L_out_extra', 0.0) or 0.0),
               outlet_stretch_extensions=coll.get('extensions'), pairs_below_floor=coll.get('pairs_below_floor'),
               part_clearance=coll['part_clearance'], collisions=coll['why'])
    if not coll['ok']:
        opts['layout_refusal'] = ('the drawn parts collide: %s.  Remedies: %s' % (
            '; '.join(coll['why']), layout_remedy(sp, sc_f, kind, n_legs, opts, coll)))
        say(opts['layout_refusal'])
    return out


def plateau_preflight(spec, opts):
    """Evaluate the kappa plateau gate before the head's sealing search (third-party review 1 of rc7, section E;
    review8 checked the equivalence on every shipped record).

    The gate's u_min is D over the widest center width, because build.gates takes u = D/w_c row by row, so its minimum
    is at the widest row. No plan draws a center wider than w_c0 (head rows draw at w_c0, and a taper only narrows), so
    u_min = D/w_c0, which equals the post-plan gate's value on all 50 drawn rc7 records. In RC7 the gate was reached
    only after the head's sealing search, which far below the plateau and with a wide chip could run for more than
    30 minutes (review 1, D 0.5 um on a 100 um center channel and a 25 x 75 mm chip). Strict mode with a head refuses
    here, giving the refusal that the post-plan gate gives, only earlier. Designs without a head and the other modes go
    on to the plan, whose envelope reports the plateau as before."""
    if opts['mode'] != 'strict' or spec.head == 'none':
        return
    u0 = S.KAPPA_G * spec.delta_star / spec.wc0
    if u0 < S.U_LO - 1e-12:
        opts['preflight_refusal'] = '; '.join([m for m in (opts.get('preflight_refusal'),
                                               'kappa plateau (u_min %.4f < %.3f at the widest centre width %.1f um, '
                                               'before the head\'s sealing search)' % (u0, S.U_LO, spec.wc0)) if m])


def _prepare(spec, opts, say):
    """Do the work of prepare() on the engine spec and its options; prepare() turns a planning failure (ValueError)
    into a SpecError."""
    notes = list(opts.get('notes_spec', []))             # what spec.to_engine found worth saying (none of it refused)
    gamma_gs_check(spec, opts, notes)
    wmin_check(spec, opts, notes)
    plateau_preflight(spec, opts)                         # the kappa plateau before the sealing search
    if spec.head != 'none' and 'kappa plateau' in (opts.get('preflight_refusal') or ''):
        # The design is refused before the head's sealing search, and the record keeps the plan of the head as given
        # (unsealed, as a seal refusal does) for its numbers.
        sc_u, _, _ = _plan(spec)
        return spec, opts, sc_u, (design_gates(spec, sc_u) if sc_u is not None else None), [], notes
    if spec.head != 'none':
        spec = seal_policy(spec, opts, notes, say)
        if not opts.get('seal_refusal'):
            spec = shrink_head(spec, opts, notes, say)
    if spec.head == 'none' and opts.get('start') == 'floor':
        spec = floor_start(spec, opts, notes, say)
    if opts.get('seal_refusal'):
        # For a refused design the record keeps the unsealed plan, if there is one, for its numbers. The envelope gates
        # of an unsealed plan are not the design's (gap 1 raised), so they are recorded without being refusal reasons.
        sc_u, _, _ = _plan(spec)
        return spec, opts, sc_u, (design_gates(spec, sc_u) if sc_u is not None else None), [], notes
    sc = S.plan(spec)
    wmin_check(spec, opts, notes, sc)
    if opts.get('F_tol') is not None and abs(float(sc.F) - float(spec.F_target)) > float(opts['F_tol']):
        # The march stops at the first row that reaches F, so the planned F can overshoot by one row's extraction.
        # When the overshoot exceeds the tolerance, the build fails Q2 F, and the note says so before drawing (review8's
        # stress test, S05).
        notes.append('F: the plan reaches F %.5f against the target %.5f (%+.5f, beyond checks.F_tolerance %g): the build '
                     'is expected to FAIL Q2 F; physics.trim_F: true trims the last gap to the target' % (
                         float(sc.F), float(spec.F_target), float(sc.F) - float(spec.F_target), float(opts['F_tol'])))
        say(notes[-1])
    problems = []                                         # spec errors, collected and raised with every other finding
    if sc.NM > 0 and sc.N < sc.NM + 2:                    # F reached inside the head (the plan-stage lattice)
        problems.append('F %.4f is reached inside the meander head (%d of the %d rows are head rows): a head needs a '
                        'linear section after it; use a smaller head (meander.head_width_um), no head (meander.enabled: '
                        'false), or a larger F' % (spec.F_target, sc.NM, sc.N))
    elif sc.NM > 0:
        try:                                              # every head unit must be drawable (smallest loop)
            L.head_top_x(sc)
        except (ValueError, IndexError) as e:
            problems.append('meander head: %s' % (e if isinstance(e, ValueError) else
                                                  'its last unit cannot be joined to the linear section (%s)'
                                                  % type(e).__name__))
    if spec.head == 'none':
        ws1 = float(np.nanmin(sc.ws[:sc.N - 1]))
        if ws1 < spec.F_min - 1e-9:
            problems.append('a sealed start without a meander head draws a %.2f um first side channel, below the %.1f um '
                            'process floor: enable a meander head (meander.enabled: true), use start.style: floor '
                            "(Gifford's rule; it over-extracts at the entrance, rows declared only in gates.mode: "
                            "reproduce: a new design usually FAILS Q2 on them), or use a finer process"
                            % (ws1, spec.F_min))
    if not problems:
        outlet_precheck(spec, sc, opts, notes, say)
        if opts.get('layout_problem'):
            problems.append(opts['layout_problem'])
    if spec.R_ext_c > 0 or spec.R_ext_s > 0:
        notes.append('outlet tubing enters the outlet closure as design values (centre %.3g, side %.3g Pa s/m^3 at water\'s '
                     'viscosity, as the chip\'s own resistances, so their ratio holds for any Newtonian fluid): it is '
                     'not on the mask and is not measured; the record reports how the last rows depend on it'
                     % (spec.R_ext_c, spec.R_ext_s))
    notes.append('measured on the drawing: every side segment, head unit length and width, U-turn side and centre arc, '
                 'outlet lead, gap opening and centre width.  Taken from the design (not measured): the post length across '
                 'the flow (gap length l_g) and so the side channels\' inner edge, the radii of the head and outlet arcs '
                 '(their widths are measured), the head rows\' centre-path lengths, and any off-chip tubing.')
    g0 = design_gates(spec, sc)
    env_bad = [k for k, ok in (('kappa plateau (u_min %.4f < %.3f)' % (g0['plateau_u_min'], S.U_LO), g0['plateau_ok']),
                               ('gap dominance M13 %.4f >= 0.05' % g0['M13_gap_dominance_max'], g0['M13_ok']),
                               ('Reynolds %.2f > 10' % g0['Re_max'], g0['Re_ok']),
                               ('Stokes gates (Re_gap %.2e at row %d, Re_delta %.2e at row %d; both must be <= %g)'
                                % (g0['Re_gap_max'], g0['Re_gap_row'], g0['Re_delta_max'], g0['Re_delta_row'], BU.RE_STOKES),
                                g0['stokes_ok']),
                               ('below the foot of the kappa curve (u_min %.4f < %.3f): no cutoff band either'
                                % (g0['plateau_u_min'], S.U_FF), g0['plateau_u_min'] >= S.U_FF),
                               ('thin lamina E11 (delta/h %.4f on the cutoff rows, %.4f on every planned row; > %g)'
                                % (g0['lamina_over_depth'], g0['lamina_over_depth_planned_max'], BU.EPS_LAM),
                                g0['lamina_ok'] and g0.get('lamina_all_rows_ok', True)),
                               ('taper rate M16 (|w_c(i) - w_c(i-1)|/(G + g_p) %.3e at row %d > %g)'
                                % (g0['taper_rate_max'], g0['taper_rate_row'], BU.EPS_TAPER), g0['taper_ok'])) if not ok]
    if problems:                                          # every finding before drawing, rather than the first one
        also = ([opts['preflight_refusal']] if opts.get('preflight_refusal') else []) + \
               ([opts['outlet_refusal']] if opts.get('outlet_refusal') else []) + \
               (['outside the validity envelope: ' + '; '.join(env_bad)] if env_bad else [])
        raise SP.SpecError('; '.join(problems) + (('.  Also found before drawing: ' + '; '.join(also)) if also else ''))
    return spec, opts, sc, g0, env_bad, notes


def _write_refusal(d, name, outdir, spec, sc, g0, notes, opts, env_bad, say, extra=None):
    """Write the refusal record (exit 5) and return it. run() writes it for a refusal found before drawing, for an
    outlet that turns out not to close while the layout is fitted or drawn, and for an outline whose parts collide,
    which is found on the outline in memory before the photomask is written (`extra` carries its collision checks)."""
    strict_env = bool(env_bad) and opts['mode'] == 'strict'
    refused, reasons, tag = refusal(opts, env_bad)
    advice = []
    if opts.get('layout_refusal'):
        advice.append('the parts of the drawn outline (the inlet with the head, each leg, each U-turn, the outlet) must '
                      'keep the process floor from every part they are not joined to, with the outlet ports separated and '
                      'at least 2 mm from the inlet port; no photomask is written for an outline that does not (the '
                      'remedies are named above).')
    if opts.get('outlet_refusal'):
        advice.append('the outlet\'s lead closes the last row (its resistance sets the last gaps\' flow): the requested '
                      'outlet cannot give the resistance or the clearances named above at this F and layout; use '
                      'outlet.style teardrop or separate, another layout or leg count, or a lower F.')
    if opts.get('seal_refusal'):
        advice.append('meander.seal: require asks for gap 1 at f*; the requested head cannot give it here (no other head '
                      'is substituted).')
    if opts.get('preflight_refusal') and 'gamma_GS' in opts['preflight_refusal']:
        advice.append('strict mode refuses a meander head whose fold width over gap lies outside gamma_GS (spec v5 '
                      'sec. 14.4, 18; CIF_v2 design()): set meander.channel_um inside [%g, %g] x process.gap_um, or build '
                      'with gates.mode: evaluate to draw it and report the departure.' % BU.GAMMA_GS)
    if opts.get('preflight_refusal') and 'w_min' in opts['preflight_refusal']:
        advice.append('strict mode refuses a gap, meander channel, centre channel or straight side channel narrower than CIF_v2\'s '
                      'w_min = max(F_min, h/AR_max) (spec v5 sec. 18): widen it, use a shallower process, or set '
                      'process.aspect_ratio_max to the deepest channel your process makes; or build with gates.mode: '
                      'evaluate to draw it and report the departure.')
    if opts.get('preflight_refusal') and 'kappa plateau' in opts['preflight_refusal']:
        # The plateau pre-flight's refusal carries the plateau's advice (review8's patch review, m1).
        advice.append('strict mode refuses a design below the kappa plateau (u = D/w_c >= %g, App. B thesis scope), here '
                      'before the head\'s sealing search: at D %.3g um the centre channel must be at most %.1f um wide '
                      '(D/%g); or build with gates.mode: evaluate to draw it and report the departure.'
                      % (S.U_LO, S.KAPPA_G * spec.delta_star, S.KAPPA_G * spec.delta_star / S.U_LO, S.U_LO))
    if strict_env:
        advice.append('strict mode refuses designs outside the validity envelope (App. B thesis scope: the kappa '
                      'plateau and the M14 channel and Stokes gates; with them the M13 gap-dominance gate of spec v5 '
                      'sec. 19). %s. Build it with gates.mode: evaluate to draw and measure the mask.'
                      % (('The cutoff is still predicted: ' + band_txt(g0['cutoff_band_um']))
                         if g0['cutoff_band_um'].get('band_um') else
                         'No cutoff is predicted there: ' + g0['cutoff_band_um'].get('basis', '')))
    rec = dict(design=dict(name=name, N=int(sc.N) if sc is not None else None, NM=int(sc.NM) if sc is not None else None,
                           F=float(sc.F) if sc is not None else None, D_star_um=S.KAPPA_G * spec.delta_star,
                           kernel_pair=BU.KERNEL_PAIR),
               refused=dict(reasons=reasons, gates=g0, advice=' '.join(advice),
                            gates_note=('the gates above are those of the unsealed plan (gap 1 raised), not of a '
                                        'design that could be built' if opts.get('seal_refusal') else
                                        'the gates above are those of the head as given, planned without its sealing '
                                        'search (the kappa plateau pre-flight refused first)'
                                        if 'kappa plateau' in (opts.get('preflight_refusal') or '') else None)),
               assessment=dict(verdict='REFUSED (%s): %s' % (tag, '; '.join(reasons)), exit_code=5),
               notes=notes, input=d, code=CODE, environment=ENV)
    if extra:
        rec['refused'].update(extra)
    json.dump(rec, open(os.path.join(outdir, name + '.json'), 'w'), indent=1, default=float)
    say(rec['assessment']['verdict'])
    return rec


def name_mask(outdir, stem, exit_code):
    """Name a photomask after its verdict and return the file name (the clearance guarantee).

    The build writes the photomask as <stem>_UNASSESSED.dxf, so that a run that stops before its verdict leaves no
    deliverable name (code review 4, F4). A build that fails (exit 2) keeps the mask for diagnosis as
    <stem>_FAILED.dxf, and never under the deliverable name <stem>.dxf, whereas any other verdict names it <stem>.dxf
    (a re-assessment that clears a FAIL renames it back). The file just assessed (<stem>_UNASSESSED.dxf) is always the
    one named, and a file of the other verdict left from an earlier run is removed."""
    ok_p, bad_p = os.path.join(outdir, stem + '.dxf'), os.path.join(outdir, stem + '_FAILED.dxf')
    un_p = os.path.join(outdir, stem + '_UNASSESSED.dxf')
    if os.path.exists(un_p):
        # The file just assessed is the one named (third-party reviews of rc7; rc7 LIMITATIONS E). RC7 left an earlier
        # <stem>.dxf in place beside it (a pass) and moved an earlier <stem>.dxf over it (a FAIL), a case that was
        # reachable through compose.array into a folder holding an earlier array of the same name.
        target, other = (bad_p, ok_p) if exit_code == 2 else (ok_p, bad_p)
        if os.path.exists(other):
            os.remove(other)
        os.replace(un_p, target)
        return os.path.basename(target)
    if exit_code == 2:                                   # no fresh file; a re-assessment renames the mask it assessed
        if os.path.exists(ok_p):
            os.replace(ok_p, bad_p)
        return stem + '_FAILED.dxf'
    if os.path.exists(bad_p) and not os.path.exists(ok_p):
        os.replace(bad_p, ok_p)
    return stem + '.dxf'


def run(user, outdir=None, verbose=True, render=True):
    """Build one device from a spec into outdir and return its record (the module docstring lists the steps)."""
    t0 = time.time()
    S.clear_caches()                                     # no state carried over from an earlier build (BUILD-1)
    d = SP.load(user)
    name = d.get('name', 'device')
    outdir = outdir or os.path.join(RUNS, name)
    os.makedirs(outdir, exist_ok=True)
    import re
    own = re.compile(re.escape(name) + r'(\.(json|dxf|pkl)|_(summary\.md|layout\.png|zooms\.png|physics\.png|audit\.json'
                     r'|FAILED\.dxf|UNASSESSED\.dxf|refused_layout\.png)'
                     r'|_x\d+(\.(json|dxf|png)|_[a-z]+\.(json|png)|_FAILED\.dxf|_UNASSESSED\.dxf))$')
    for f in os.listdir(outdir):                         # a record must never mix runs (this name's files only)
        if own.match(f):
            os.remove(os.path.join(outdir, f))

    def say(m):
        if verbose:
            print('[%s] %s' % (name, m), flush=True)
    spec, opts, sc, g0, env_bad, notes = prepare(d, say)
    refused, reasons, tag = refusal(opts, env_bad)
    if refused:
        return _write_refusal(d, name, outdir, spec, sc, g0, notes, opts, env_bad, say)
    if env_bad:
        notes.append('outside the validity envelope: %s (built because gates.mode = %s)' % ('; '.join(env_bad), opts['mode']))
    try:
        kind, n_legs, spec, sc = fit_layout(spec, sc, opts, notes, say)
        rep, fluid, sc_f, info = BU.build(spec, kind, n_legs, outdir, name, failsafe=opts['failsafe'],
                                          do_audit=opts['drc'], verbose=verbose, R_turn=opts['leg_spacing'],
                                          cut_tol=opts['cut_tol'], mode=opts['mode'], seal=opts['seal'],
                                          ar_max=opts.get('AR_max', BU.AR_MAX_DEFAULT))
        if getattr(spec, 'L_out_extra', 0.0):            # the outlet stretch was lengthened (the clearance guarantee)
            notes.append('layout: the outlet met the inlet end, so its stretch was lengthened by %.0f um (the closure '
                         'takes the longer stretch; build.plan_and_draw)' % spec.L_out_extra)
    except L.OutletOutsideTables as e:                    # the stretch-end tables' domain (strict mode)
        opts['outlet_refusal'] = ('the %s outlet: %s: build it in evaluate or reproduce mode, where the end values '
                                  'are held and reported, or change the outlet' % (opts.get('outlet'), e))
        return _write_refusal(d, name, outdir, spec, sc, g0, notes, opts, env_bad, say)
    except L.OutletInfeasible as e:                       # found only while fitting or drawing, with the same refusal
        opts['outlet_refusal'] = 'the %s outlet cannot close this design (F %.4f): %s' % (opts.get('outlet'),
                                                                                          spec.F_target, e)
        return _write_refusal(d, name, outdir, spec, sc, g0, notes, opts, env_bad, say)
    except BU.LayoutCollision as e:                       # refused before anything is written (the clearance guarantee)
        opts['layout_refusal'] = ('the drawn parts collide: %s.  Remedies: %s' % (
            '; '.join(e.coll['why']), layout_remedy(spec, e.sc, kind, n_legs, opts, e.coll)))
        rec = _write_refusal(d, name, outdir, spec, sc, g0, notes, opts, env_bad, say,
                             extra=dict(layout=dict(layout=kind, legs=n_legs, collisions=e.coll['why'],
                                                    part_clearance=e.coll['part_clearance'],
                                                    pairs_below_floor=e.coll.get('pairs_below_floor'),
                                                    port_spacing=e.coll['port_spacing'], separation=e.coll['separation'],
                                                    outlet_stretch_extensions=e.coll.get('extensions'),
                                                    holes=len(e.coll['holes']), holes_expected=int(e.coll['exp_holes']))))
        if render:
            try:
                rec['files'] = dict(figures=[RN.refused_layout(e.fluid, e.info, e.coll, outdir, name,
                                                               rec['assessment']['verdict'])])
                json.dump(rec, open(os.path.join(outdir, name + '.json'), 'w'), indent=1, default=float)
            except Exception as e_:                       # a convenience figure; the record says why it is missing
                rec['files'] = dict(figures=[], figure_error='%s: %s' % (type(e_).__name__, e_))
                json.dump(rec, open(os.path.join(outdir, name + '.json'), 'w'), indent=1, default=float)
        return rec
    chip_rect = None
    if opts['chip'][0] and opts['chip'][1]:
        chip_rect = (opts['chip'][1] * 1000.0, opts['chip'][0] * 1000.0)
    rep['assessment'] = assess(rep, spec, opts, fluid, notes, chip_rect)
    rep['code'] = CODE                                   # hashed when this process imported the package
    rep['environment'] = ENV                             # the test suite and the environment (ADJ-C7)
    rep['input'] = d
    figs = []
    if render:
        import pickle
        arr = pickle.load(open(os.path.join(outdir, name + '.pkl'), 'rb'))['arrays']
        chip = None
        if chip_rect:
            b = fluid.bounds
            cx, cy = 0.5 * (b[0] + b[2]), 0.5 * (b[1] + b[3])
            chip = (chip_rect[0], chip_rect[1], (cx - chip_rect[0] / 2, cy - chip_rect[1] / 2))
        title = '%s: D %.3f µm, F %.4f; %s, %s taper, %s head, %s inlet, %s outlet; %s' % (
            name, S.KAPPA_G * spec.delta_star, spec.F_target, kind + ('' if n_legs == 1 else ' ×%d' % n_legs),
            spec.taper, (opts['style'] or 'no'), opts['inlet'], opts['outlet'], rep['assessment']['verdict'])
        rep['title'] = title
        figs = RN.all_figures(fluid, sc_f, info, arr, rep, outdir, name, title, chip)
    from opencif.geometry import dxf as io_dxf
    dxf_name = name_mask(outdir, name, rep['assessment'].get('exit_code'))
    dxf_path = os.path.join(outdir, dxf_name)
    drc_ = (rep.get('validation') or {}).get('drc')
    if isinstance(drc_, dict) and drc_.get('file') and drc_['file'] != dxf_name:
        # 'file' is the name that the DRC opened (<name>_UNASSESSED.dxf, before the verdict), and the same bytes now
        # carry the verdict's name (rc7 LIMITATIONS E; third-party review 2, C1).
        drc_['file_final'] = dxf_name
    pa_ = (rep.get('validation') or {}).get('physics_audit')
    if isinstance(pa_, dict) and (pa_.get('photomask') or {}).get('file') not in (None, dxf_name):
        pa_['photomask']['file_final'] = dxf_name         # the audited file (the same bytes) under the verdict's name
    rep['files'] = dict(dxf=dxf_name, dxf_sha256=io_dxf.sha256(dxf_path) if os.path.exists(dxf_path) else None,
                        record=name + '.json', figures=[os.path.basename(f) for f in figs])
    if d.get('array'):
        from . import compose as CP
        ar = d['array']
        json.dump(rep, open(os.path.join(outdir, name + '.json'), 'w'), indent=1, default=float)
        try:
            rep['array'] = CP.array(outdir, name, int(ar['count']), outdir, ar.get('pitch_um'), ar.get('inlet', 'shared'),
                                    ar.get('bus_width_um'), chip=chip_rect, margin_um=opts['margin'], do_drc=opts['drc'])
            say('array: %s' % rep['array']['verdict'])
        except ValueError as e:
            rep['array'] = dict(verdict='NOT BUILT: %s' % e)
            say('array not built: %s' % e)
    rep['seconds_total'] = time.time() - t0
    json.dump(rep, open(os.path.join(outdir, name + '.json'), 'w'), indent=1, default=float)
    summary_md(rep, os.path.join(outdir, name + '_summary.md'))
    say('%s in %.0f s' % (rep['assessment']['verdict'], time.time() - t0))
    return rep


def reassess(outdir, name, render=True):
    """Re-run the analysis (Q1, Q2, declared rows, verdict, summary, figures) on a finished build's saved arrays and
    record, without drawing again. It is used when the assessment rules change, and the measured networks are not
    recomputed."""
    import pickle
    rep = json.load(open(os.path.join(outdir, name + '.json')))
    pk = pickle.load(open(os.path.join(outdir, name + '.pkl'), 'rb'))
    arr, fluid = pk['arrays'], pk['fluid']
    d = rep['input']
    spec_u, opts = SP.to_engine(d)
    spec = arr['spec']
    V = rep['validation']
    if opts['seal'] == 'auto' and opts['mode'] != 'strict' and not (V.get('components_drawn') or {}).get('gap1_sealed', True):
        opts['seal'] = 'allow'                            # as seal_policy did when the record was built
    if 'stokes_ok' not in V.get('gates', {}):
        # A record from before the Stokes gates has its gates recomputed from the saved design arrays. M13 needs the
        # closure resistances, which are not saved, so its value is kept.
        import types
        f_ = np.asarray(arr['f'], float)
        P_ = np.asarray(arr['P'], float)
        wc_ = np.asarray(arr['wc'], float)
        Pin_ = np.concatenate([[1.0], P_[:-1]])
        wL_ = np.array([spec.G if k == 'L' else (spec.w_arm if k == 'A' else ((spec.w_narrow or spec.GS) if k == 'N'
                        else spec.GS)) for k in arr['kind']])
        dl_ = np.array([S.x_of(fi) * wi if 0 < fi < 0.5 else float('nan') for fi, wi in zip(f_, wc_)])
        ns = types.SimpleNamespace(wc=wc_, P=P_, f=f_, Qg=f_ * Pin_, wL=wL_, delta=dl_, Rg=np.zeros_like(f_),
                                   Rc=np.ones_like(f_))
        g_new = BU.gates(spec, ns, S.KAPPA_G * spec.delta_star)
        for k_ in ('M13_gap_dominance_max', 'M13_ok'):
            g_new[k_] = V['gates'][k_]
        V['gates'] = g_new
    D = S.KAPPA_G * spec.delta_star
    nw = dict(fA=np.asarray(arr['fA']), fB=np.asarray(arr['fB']), dA=np.asarray(arr['dA']), dB=np.asarray(arr['dB']),
              F=V['network']['F_as_drawn'])
    inten = dict(fA=np.asarray(arr['fA_int']), fB=np.asarray(arr['fB_int']), dA=np.asarray(arr['dA_int']),
                 dB=np.asarray(arr['dB_int']), F=V['network']['F_intended'])
    fsc = np.asarray(arr['fstar_cut'])
    V['drawing_fidelity'] = BU.drawing_fidelity(nw, inten, fsc, int(arr['N']), D)
    nets = {k: dict((kk, np.asarray(vv) if kk != 'F' else vv) for kk, vv in v.items())
            for k, v in (arr.get('nets') or dict(full=inten)).items()}
    decl = BU.declared_rows(spec, np.asarray(arr['f']), fsc, list(arr['turns']), D, opts['cut_tol'], nets, opts['mode'],
                            opts['seal'])
    V['realised_cutoff'] = BU.realised_cutoff(None, nw, D, opts['cut_tol'], decl)
    V['intended_cutoff'] = BU.realised_cutoff(None, inten, D, opts['cut_tol'], decl)
    arr['declared'] = decl
    chip_rect = (opts['chip'][1] * 1000.0, opts['chip'][0] * 1000.0) if (opts['chip'][0] and opts['chip'][1]) else None
    notes = rep.get('assessment', {}).get('notes', [])
    rep['assessment'] = assess(rep, spec, opts, fluid, notes, chip_rect)
    rep['reassessed'] = dict(code=code_stamp(), note='checks recomputed from the saved arrays; the drawing, its measurement '
                             'and the intended network are those of the build stamped in `code`')
    rep['title'] = '; '.join((rep.get('title') or name).split('; ')[:2]) + '; ' + rep['assessment']['verdict']
    if rep.get('files'):                                  # a FAIL's mask is <name>_FAILED.dxf (the clearance guarantee)
        from opencif.geometry import dxf as io_dxf
        fn_ = name_mask(outdir, name, rep['assessment'].get('exit_code'))
        rep['files']['dxf'] = fn_
        drc_ = (rep.get('validation') or {}).get('drc')
        if isinstance(drc_, dict) and drc_.get('file') and drc_['file'] != fn_:
            drc_['file_final'] = fn_                      # as in run()
        pa_ = (rep.get('validation') or {}).get('physics_audit')
        if isinstance(pa_, dict) and (pa_.get('photomask') or {}).get('file') not in (None, fn_):
            pa_['photomask']['file_final'] = fn_          # as in run()
        rep['files']['dxf_sha256'] = io_dxf.sha256(os.path.join(outdir, fn_)) if os.path.exists(os.path.join(outdir, fn_)) else None
    with open(os.path.join(outdir, name + '.pkl'), 'wb') as fh:
        pickle.dump(dict(fluid=fluid, arrays=arr), fh)
    json.dump(rep, open(os.path.join(outdir, name + '.json'), 'w'), indent=1, default=float)
    summary_md(rep, os.path.join(outdir, name + '_summary.md'))
    if render:
        RN.rerender(outdir, name)
    return rep
