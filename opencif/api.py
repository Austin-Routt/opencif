"""opencif.api: one call from the two decoupled targets to a validated photomask and its figures.

    import opencif
    r = opencif.generate(4.0, 0.5)                                   # D = 4 um, F = 0.5, every other choice defaulted
    r = opencif.generate(2.5, 0.6, process='dinh2024', chip=(25, 75), layout='serpentine',
                         taper=('constant_shear', 70, 50), meander='pack', inlet='filter', outlet='teardrop')
    r.verdict, r.exit_code, r.dxf, r.figures, r.record                 # r.ok is True for exit code 0
    r.over_chip_area, r.over_aspect_ratio      # the chip area is reported and is never part of r.ok, whereas
                                               # the aspect ratio is a pass condition in strict mode and is
                                               # reported in the other modes

The size cutoff D (cutoff_um) and the filtrate fraction F (or the filtrate-to-retentate flow ratio F/(1-F)) are the
only two design targets. The depth, the gap, the post sizes, the flow, and the fluid (water by default) are model
inputs, which set the mask and the validity checks rather than the targets. Every other argument chooses a part or the
fabrication context, and each maps onto exactly one section of the YAML spec, as the table below shows
(docs/OPTIONS.md lists them all, with their codes and defaults):

    argument        spec section        accepts
    process         process             a name in PROCESSES or a dict (depth_um, gap_um, min_feature_um, ...)
    chip            chip                (width_mm, length_mm[, margin_um]) or a dict; None = no area
    layout, legs    layout              'linear' | 'serpentine' | 'paperclip' (or a dict); legs: int or 'auto'
    taper           centre_channel      'constant' | 'constant_shear' | (style, width_in_um[, width_min_um]) | a dict
    meander         meander             None/False (no head) | 'pack' | 'arms' | 'hairpins' | 'narrow' | a dict
    seal            meander.seal        'require' | 'auto' | 'allow'
    inlet           inlet               'plain' | 'filter' | 'dinh_filter' (or a dict)
    outlet          outlet              'separate' | 'teardrop' | 'ring' (or a dict)
    start           start               None | 'sealed' | 'floor' (or a dict; no-head devices only)
    flow_uL_min     operation           the operating flow (Reynolds and Stokes gates)
    mode            gates.mode          'strict' | 'evaluate' | 'reproduce'
    array           array               None | an int (count) | a dict (count, pitch_um, inlet, bus_width_um)
    overrides       any section         a dict merged last, section by section (physics, checks, ...)

The call writes the spec it built to <out>/<name>_spec.yaml, so the same device can be rebuilt from the command line
with `python -m opencif build <out>/<name>_spec.yaml`. Nothing here adds physics, because the call builds the spec
dict and hands it to pipeline.run, whose schema check, gates, and validators apply unchanged.
"""
import copy
import os
from dataclasses import dataclass, field

# The named processes. Each is the process block of a spec in this repository (`PROCESS_SOURCES` gives the source),
# and a dict can always be passed instead. The keys are the depth h, the gap G, the post length g_p (along the flow),
# the post width l_g (across the flow), and the process floor.
PROCESSES = {
    'dinh2024': dict(depth_um=140.0, gap_um=19.0, post_length_um=26.1, post_width_um=54.0, min_feature_um=19.0),
    'xia2016': dict(depth_um=140.0, gap_um=16.0, post_length_um=35.8, post_width_um=70.9, min_feature_um=16.0,
                    min_radius_um=4.5),
    'gifford2014': dict(depth_um=80.0, gap_um=20.5, post_length_um=25.3, post_width_um=59.5, min_feature_um=19.0),
    'generic19': dict(depth_um=100.0, gap_um=19.0, post_length_um=27.0, post_width_um=50.0, min_feature_um=19.0),
}
PROCESS_SOURCES = {
    'dinh2024': 'presets/dinh2024_v1_specified.yaml (Dinh 2024 p. 917, Fig. 2 caption: 19 um gaps, 140 um deep)',
    'xia2016': 'presets/xia2016_rs7_specified.yaml (Xia 2016 RS-7 as measured on its figures)',
    'gifford2014': 'presets/gifford2014_fig4_design.yaml with a 19 um floor, as examples/gifford_targets_meander_head_19um.yaml',
    'generic19': "examples/matrix/*.yaml (this project's 19 um, 100 um deep process)",
}
DEFAULT_U = 0.05            # the default center width is w_c = D/u, and u = 0.05 sits on the kappa plateau (u >= 0.035)


@dataclass
class Result:
    """The result of generate(), which holds the record and the paths of its files."""
    record: dict
    spec: dict
    out: str
    verdict: str = ''
    exit_code: int = None
    dxf: str = None
    record_path: str = None
    summary: str = None
    spec_path: str = None
    figures: list = field(default_factory=list)

    @property
    def ok(self):
        """True for exit code 0, which means that every hard check of the generator, the physics audit of the photomask
        file among them, passed inside the validity envelope. This is a statement about the 1-D model on the drawn
        mask, and it does not say that the chip is fabricable or that a device was made or measured. The chip area is
        reported and is never a pass condition (App. B, claim C8, CQ5), whereas the aspect ratio is a pass condition in
        strict mode and is reported in the other modes (docs/LIMITATIONS.md, item 7). See over_chip_area,
        over_aspect_ratio, and record['assessment']['status']."""
        return self.exit_code == 0

    @property
    def physics_audit(self):
        """The status of the physics audit of the photomask file (opencif.checks.auditor), which is 'PASS',
        'PASS_OUTSIDE_ENVELOPE', or 'FAIL', or None for a refused design, which has no photomask. Every check of the
        audit is in record['validation']['physics_audit']."""
        return (((self.record or {}).get('validation') or {}).get('physics_audit') or {}).get('status')

    @property
    def over_chip_area(self):
        """True when the drawn device exceeds the chip area given in the spec, which is reported and stated in the
        verdict."""
        return bool(((self.record or {}).get('assessment') or {}).get('over_chip_area', False))

    @property
    def over_aspect_ratio(self):
        """True when the narrowest CIF channel is deeper than process.aspect_ratio_max allows. This is a FAIL in strict
        mode, whereas evaluate and reproduce modes report it and say so in the verdict."""
        return bool(((self.record or {}).get('assessment') or {}).get('over_aspect_ratio', False))


def _section(value, key, name):
    """Turn a part argument, which is a style string or a dict for the whole section, into its spec section."""
    if value is None:
        return None
    if isinstance(value, dict):
        return copy.deepcopy(value)
    if isinstance(value, str):
        return {key: value}
    raise TypeError('%s: a style string or a dict, got %r' % (name, value))


def _num(x):
    """Format a number for a device name, with 'p' for the decimal point and 'm' for the minus sign."""
    return ('%g' % x).replace('.', 'p').replace('-', 'm')


def build_spec(cutoff_um, filtrate_fraction=None, *, filtrate_to_retentate_ratio=None, process='generic19', chip=None,
               layout='linear', legs=None, taper='constant', width_in_um=None, width_min_um=None, meander='pack',
               seal=None, inlet='plain', outlet='separate', start=None, flow_uL_min=None, mode='strict', array=None,
               name=None, description=None, overrides=None):
    """Return the spec dict that generate() builds. No file is written, and nothing is validated beyond what
    spec.load checks."""
    from .spec import SpecError                         # one exception type for a bad input (rc7 raised ValueError)

    def num(v, field):
        # A flag is refused as a number, NumPy's bool included (RC8 code review, B1, N1).
        if isinstance(v, bool) or getattr(getattr(v, 'dtype', None), 'kind', '') == 'b':
            raise SpecError('%s must be a number, not %r' % (field, v))
        try:
            return float(v)
        except (TypeError, ValueError):
            raise SpecError('%s must be a number, got %r' % (field, v))
    D = num(cutoff_um, 'targets.cutoff_um')
    if filtrate_fraction is not None:
        filtrate_fraction = num(filtrate_fraction, 'targets.filtrate_fraction')
    if filtrate_to_retentate_ratio is not None:
        filtrate_to_retentate_ratio = num(filtrate_to_retentate_ratio, 'targets.filtrate_to_retentate_ratio')
    Fv = (filtrate_fraction if filtrate_fraction is not None else
          (filtrate_to_retentate_ratio / (1.0 + filtrate_to_retentate_ratio) if filtrate_to_retentate_ratio else 0.0))
    d = dict(name=name or 'cif_D%s_F%s' % (_num(D), _num(round(Fv, 4))))
    t = dict(cutoff_um=D)
    if filtrate_fraction is not None:
        t['filtrate_fraction'] = float(filtrate_fraction)
    if filtrate_to_retentate_ratio is not None:
        t['filtrate_to_retentate_ratio'] = float(filtrate_to_retentate_ratio)
    d['targets'] = t
    if isinstance(process, str):
        if process not in PROCESSES:
            raise SpecError('process %r: one of %s, or a dict' % (process, ', '.join(PROCESSES)))
        d['process'] = dict(PROCESSES[process])
    else:
        d['process'] = copy.deepcopy(process)
    if chip is not None:
        if isinstance(chip, dict):
            d['chip'] = copy.deepcopy(chip)
        else:
            c = list(chip)
            d['chip'] = dict(width_mm=num(c[0], 'chip.width_mm'), length_mm=num(c[1], 'chip.length_mm'))
            if len(c) > 2:
                d['chip']['margin_um'] = num(c[2], 'chip.margin_um')
    lay = _section(layout, 'style', 'layout')
    if legs is not None:
        lay['legs'] = legs
    d['layout'] = lay
    notes = []
    if isinstance(taper, (tuple, list)):
        cc = dict(taper=taper[0], width_in_um=num(taper[1], 'centre_channel.width_in_um'))
        if len(taper) > 2:
            cc['width_min_um'] = num(taper[2], 'centre_channel.width_min_um')
    else:
        cc = _section(taper, 'taper', 'taper')
    if width_in_um is not None:
        cc['width_in_um'] = num(width_in_um, 'centre_channel.width_in_um')
    if width_min_um is not None:
        cc['width_min_um'] = num(width_min_um, 'centre_channel.width_min_um')
    if 'width_in_um' not in cc:
        cc['width_in_um'] = round(D / DEFAULT_U, 1)
        notes.append('centre_channel.width_in_um = D/%.2f = %.1f um (generate default: u = D/w_c = %.2f, on the kappa '
                     'plateau)' % (DEFAULT_U, cc['width_in_um'], DEFAULT_U))
    if cc.get('taper') == 'constant_shear' and 'width_min_um' not in cc:
        cc['width_min_um'] = round(0.7 * cc['width_in_um'], 1)
        notes.append('centre_channel.width_min_um = 0.7 x width_in_um = %.1f um (generate default: Dinh 2024 narrows '
                     '75 -> 52 um, 0.69)' % cc['width_min_um'])
    d['centre_channel'] = cc
    if meander is None or meander is False:
        d['meander'] = dict(enabled=False)
    else:
        m = _section(meander, 'style', 'meander')
        m.setdefault('enabled', True)
        if seal is not None:
            m['seal'] = seal
        d['meander'] = m
    if seal is not None and not d['meander'].get('enabled'):
        raise SpecError('seal applies to a meander head (meander=None draws none)')
    if inlet is not None:
        d['inlet'] = _section(inlet, 'style', 'inlet')
    if outlet is not None:
        d['outlet'] = _section(outlet, 'style', 'outlet')
    if start is not None:
        d['start'] = _section(start, 'style', 'start')
    if flow_uL_min is not None:
        d['operation'] = dict(flow_uL_min=num(flow_uL_min, 'operation.flow_uL_min'))
    d['gates'] = dict(mode=mode)
    if array is not None:
        if isinstance(array, bool):                      # True is not a count (review 3, E9)
            raise SpecError('array must be a copy count (an int) or a mapping, not %r' % (array,))
        d['array'] = dict(count=int(array)) if isinstance(array, int) else copy.deepcopy(array)
    for sec, v in (overrides or {}).items():
        if isinstance(v, dict) and isinstance(d.get(sec), dict):
            d[sec].update(copy.deepcopy(v))
        else:
            d[sec] = copy.deepcopy(v)
    desc = [description] if description else []
    desc += ['made by opencif.generate()'] + notes
    d['description'] = '; '.join(desc)
    return d


def generate(cutoff_um, filtrate_fraction=None, *, out=None, render=True, drc=True, check_only=False, verbose=False,
             **parts):
    """Build a validated DXF photomask, its figures, and its record from D (in um) and F, and return them as a Result.
    The module docstring lists the part arguments (process, chip, layout, legs, taper, width_in_um, width_min_um,
    meander, seal, inlet, outlet, start, flow_uL_min, mode, array, name, description, filtrate_to_retentate_ratio, and
    overrides).

    With check_only=True, the call plans the design and evaluates its gates without writing a photomask, as
    `python -m opencif check` does. The returned record then holds N, the planned F, the seal decision, the
    feasibility of the outlet (its plan on the layout's leg counts, from pipeline.outlet_precheck), the gates, and the
    layout with the collision checks of its outline. The outline is drawn in memory exactly as a build would draw it,
    so parts that collide are refused here as a build refuses them. The exit code is the one a build would give (5 for
    a refused design, 3 for one outside the envelope, and 0 for one inside it), and the measurement of the drawn mask,
    the network, and the DRC are not run. With check_only='plan', the call plans and gates only, without the outline.
    This takes seconds for a design without a head, although parts that collide are then found only by a build or a
    full check. The record of a build is the same one that `python -m opencif build` writes, and its 'assessment'
    carries the verdict and every check."""
    from . import pipeline as PL
    from . import spec as SP
    d = build_spec(cutoff_um, filtrate_fraction, **parts)
    if not drc:
        d.setdefault('checks', {})['drc'] = False
    SP.load(d)                                           # schema and value checks before anything is planned
    name = d['name']
    out = os.path.abspath(out or os.path.join(PL.RUNS, name))
    os.makedirs(out, exist_ok=True)

    # The spec is written after the run. run() clears this name's earlier output files from outdir, but its pattern
    # (in pipeline.run) spares <name>_spec.yaml and <name>_spec.json, so the order is a precaution rather than a
    # requirement.
    def write_spec():
        try:
            import yaml
            path = os.path.join(out, name + '_spec.yaml')
            with open(path, 'w') as fh:
                yaml.safe_dump(d, fh, sort_keys=False)
        except ImportError:                              # PyYAML is optional; JSON is also a valid spec
            import json
            path = os.path.join(out, name + '_spec.json')
            json.dump(d, open(path, 'w'), indent=1)
        return path
    if check_only:
        from opencif.physics import schedule as S
        from . import build as BU
        S.clear_caches()                                  # as in pipeline.run (BUILD-1, R16)
        say_ = print if verbose else (lambda m: None)
        # The plan starts from the loaded spec, as in run() (code review 4, F9).
        spec, opts, sc, g, env_bad, notes = PL.prepare(SP.load(d), say_)
        # The layout and the outline that a build would draw are made in memory, with their collision checks (the
        # clearance guarantee).
        lay = (PL.layout_check(spec, sc, opts, notes, say_) if check_only != 'plan' and not PL.refusal(opts, env_bad)[0]
               else None)
        refused, reasons, _tag = PL.refusal(opts, env_bad)
        rec = dict(design=dict(name=name, N=int(sc.N) if sc is not None else None,
                               NM=int(sc.NM) if sc is not None else None, F_planned=float(sc.F) if sc is not None else None,
                               D_um=S.KAPPA_G * spec.delta_star, sealed=PL.is_sealed(sc) if sc is not None else False,
                               kernel_pair=BU.KERNEL_PAIR),
                   gates=g, seal_refusal=opts.get('seal_refusal'), outlet_refusal=opts.get('outlet_refusal'),
                   layout_refusal=opts.get('layout_refusal'), layout=lay,
                   preflight_refusal=opts.get('preflight_refusal'),
                   envelope_failures=env_bad, notes=notes, input=d,
                   code=PL.CODE, environment=PL.ENV)
        # A build of the planned design would give this exit code (5 refused, 3 outside the envelope, 0 inside).
        code = 5 if refused else (3 if env_bad else 0)
        state = ('planned only (no outline: check_only=\'plan\')' if check_only == 'plan' else
                 'planned and outlined in memory (no photomask written)')
        verdict = ('WOULD BE REFUSED: %s' % '; '.join(reasons) if refused else
                   '%s: %s' % (state, ('outside the envelope: ' + '; '.join(env_bad)) if env_bad else 'inside the envelope'))
        return Result(record=rec, spec=d, out=out, spec_path=write_spec(), exit_code=code, verdict=verdict)
    rec = PL.run(d, out, verbose=verbose, render=render)
    a = rec.get('assessment', {})
    files = rec.get('files') or {}
    j = lambda f: os.path.join(out, f) if f else None
    summ = os.path.join(out, name + '_summary.md')
    return Result(record=rec, spec=d, out=out, verdict=a.get('verdict', ''), exit_code=a.get('exit_code'),
                  dxf=j(files.get('dxf')), record_path=os.path.join(out, name + '.json'),
                  summary=summ if os.path.exists(summ) else None, spec_path=write_spec(),
                  figures=[j(f) for f in files.get('figures', [])])
