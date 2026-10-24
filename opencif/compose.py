"""opencif.compose: parallel arrays of a validated element.

Flow in these devices is laminar, and in the 1-D network every resistance is independent of the flow. App. A shows that
in a fully developed unidirectional duct the Navier-Stokes equation reduces exactly to the Stokes equation "whatever the
Reynolds number", and the glossary describes the flow as creeping at the gap and lamina scales, "not globally Stokes".
Thus every gap fraction f_i, and hence the cutoff D and the filtrate fraction F, depends only on the element's own
resistances and on its outlet boundary condition rather than on how much flow the element receives. An array of
identical elements whose outlets each end in their own (punched) ports at the same pressure therefore realizes the
element's validated D and F in every element, whatever the inlet manifold does, and the manifold only sets each
element's throughput, which is computed here. Merging outlets on the chip would change each element's outlet closure,
so it is not offered.

An array stacks n copies of the element across the chip at a pitch that keeps the process floor (and a clearance)
between them. With `inlet: shared`, the elements' inlet ports are joined by straight stubs to one bus channel with a
single inlet port at its middle. That requires every inlet port to be reachable from outside the element in one layer,
as in a linear or serpentine element, whereas a paperclip's inlet sits inside its loops. Otherwise each element keeps
its own inlet port.

The element must be a clean PASS whose photomask passed the physics auditor (its record's validation.physics_audit;
opencif.checks.auditor), because the copies are its geometry, translated. The checks on the drawn array are (i) the DRC
at the floor, (ii) one fluid component (shared inlet) or n (separate), (iii) in every element, a product port that is
cut off from the filtrate port(s) and the inlet once its outlet stretch is cut, (iv) a bus and stubs that touch each
element only at its own inlet port, and (v) the footprint against the chip, which is reported rather than a pass
condition (App. B C8 CQ5).
"""
import json
import math
import os
import pickle
import subprocess
import sys

import numpy as np
from shapely import affinity
from shapely.geometry import Point, box
from shapely.ops import unary_union

from . import build as BU
from opencif.physics import schedule as S  # noqa: E402
from opencif.geometry import layout as L  # noqa: E402
from opencif.geometry import dxf as io_dxf  # noqa: E402


def element_ok(rec):
    """Return (ok, reason), which says whether the element may be arrayed, because only a clean pass composes. A clean
    pass has exit code 0 and a PASS verdict that
    names no departure, so every hard check passed, no row was exempted (a reproduction's declared departures are a
    statement about a published device rather than a validated element), and the design is inside the validity
    envelope. Its photomask must also have passed the physics auditor, and it must have been built with the stokes
    bend model."""
    a = rec.get('assessment', {})
    verdict, code = a.get('verdict', ''), a.get('exit_code')
    if code != 0 or not verdict.startswith('PASS') or 'OUTSIDE' in verdict or 'departure' in verdict:
        return False, ('only a clean PASS element (exit code 0, inside the envelope, no exempted rows) can be arrayed; '
                       'this element has exit code %r and verdict %r' % (code, verdict))
    pa = (rec.get('validation') or {}).get('physics_audit') or {}
    if pa.get('status') != 'PASS':
        return False, ('only an element whose photomask passed the physics auditor can be arrayed (validation.physics_audit '
                       'status %r)' % pa.get('status'))
    bm = (rec.get('design', {}).get('spec', {}) or {}).get('bend_model', 'stokes')
    if bm != 'stokes':
        return False, 'the element was built with bend_model %r (comparison only): it cannot be arrayed' % bm
    return True, ''



def inter_copy_port_spacing(inlets, outlet_ports, outlet_radii, count, pitch, required_um):
    """Apply the element's 2 mm rule between its inlet port and every outlet port between the copies of an array (code
    review 4 F10; third-party review 2 of rc7, C4, corrected).

    Each inlet port that is an access port, which is each copy's own (separate inlets) or the shared port (`copy`
    None), is measured edge to edge against every other copy's outlet ports. The outlet ports of neighboring copies
    are reported rather than gated, because the element sets no rule between its own outlet ports. The two port
    arguments are

        inlets        [((x, y), radius, copy or None)]
        outlet_ports  {name: (x, y)} of copy 0 (copy j is shifted by j pitch in y)
    """
    gaps = [(math.hypot(float(xo) - xi, float(yo) + jo * pitch - yi) - ri - float(outlet_radii[k]), ji, jo, k)
            for (xi, yi), ri, ji in inlets for jo in range(count) if jo != ji for k, (xo, yo) in outlet_ports.items()]
    g = min(gaps) if gaps else None
    oo = [math.hypot(float(xa) - float(xb), float(ya) + pitch - float(yb)) - float(outlet_radii[ka]) - float(outlet_radii[kb])
          for ka, (xa, ya) in outlet_ports.items() for kb, (xb, yb) in outlet_ports.items()] if count > 1 else []
    return dict(ok=bool(g is None or g[0] >= required_um - 1e-6), required_um=float(required_um),
                inlet_to_other_copies_outlet_ports_min_um=(round(float(g[0]), 1) if g else None),
                closest=(dict(inlet_copy='shared' if g[1] is None else g[1], outlet_copy=g[2], port=g[3]) if g else None),
                outlet_ports_of_neighbouring_copies_min_um=(round(float(min(oo)), 1) if oo else None),
                note='edge to edge; the outlet ports of neighbouring copies are reported, not gated (the element sets no '
                     'rule between its own outlet ports): tubing for them is the user\'s to place')

def array(element_dir, name, count, outdir, pitch_um=None, inlet='shared', bus_width_um=None, clearance_um=None,
          chip=None, margin_um=500.0, do_drc=True):
    """Compose count copies of the validated element `name` from element_dir into outdir, and return the array's
    record with its photomask, checks, and verdict. A ValueError refuses an element that is not a clean pass, a pitch
    that leaves less than twice the floor, and an array whose parts merge or disconnect or whose inlet crowds another
    copy's outlet ports, before any photomask is written."""
    # This array's earlier outputs are cleared first, so that a composition into a used folder cannot leave an earlier
    # photomask, record, audit, or figure under its names. Only this array's own six names are cleared, and the element
    # and other arrays are untouched (third-party review 2 of rc7, C1; rc7 LIMITATIONS E).
    for suffix in ('.json', '.dxf', '_FAILED.dxf', '_UNASSESSED.dxf', '_layout.png', '_audit.json'):
        p_ = os.path.join(outdir, '%s_x%d%s' % (name, count, suffix))
        if os.path.isfile(p_):
            os.remove(p_)
    rec = json.load(open(os.path.join(element_dir, name + '.json')))
    pk = pickle.load(open(os.path.join(element_dir, name + '.pkl'), 'rb'))
    fl, arr = pk['fluid'], pk['arrays']
    spec = arr['spec']
    verdict = rec['assessment']['verdict']
    ok_, why = element_ok(rec)
    if not ok_:
        raise ValueError(why)
    b = fl.bounds
    H = b[3] - b[1]
    clear = clearance_um if clearance_um is not None else max(300.0, 4 * spec.F_min)
    pitch = pitch_um if pitch_um is not None else H + clear
    if pitch < H + 2 * spec.F_min:
        raise ValueError('pitch %.0f um leaves less than twice the floor between elements (element height %.0f um)' % (pitch, H))
    x_in, r_in = float(arr['inlet']['x_port_c']), float(spec.port_r_in)
    port_in = Point(x_in, 0.0).buffer(r_in, quad_segs=64)
    # A ray from the port's left edge to beyond the element tests whether the inlet port is reachable from the left in
    # one layer.
    w_bus = bus_width_um or max(4 * spec.F_min, 2 * min(r_in, 250.0))
    # The stub's own corridor (its full width plus the floor on each side) must be free of the element's fluid.
    hw_ray = w_bus / 2.0 + spec.F_min
    ray = box(b[0] - 10.0, -hw_ray, x_in - r_in + 1.0, hw_ray)
    blocked = fl.difference(port_in.buffer(1.0)).intersection(ray).area > 1.0
    shared = inlet == 'shared' and not blocked
    notes = []
    if inlet == 'shared' and blocked:
        notes.append('the inlet port is enclosed by the element (e.g. a paperclip): each element keeps its own inlet port')
    els = [affinity.translate(fl, 0.0, j * pitch) for j in range(count)]
    parts = list(els)
    bus_geo = None
    if shared:
        x_bus = b[0] - max(1000.0, 3 * w_bus)
        stubs = [box(x_bus - w_bus / 2, j * pitch - w_bus / 2, x_in, j * pitch + w_bus / 2) for j in range(count)]
        y_mid = 0.5 * (count - 1) * pitch
        bus = box(x_bus - w_bus / 2, -w_bus / 2, x_bus + w_bus / 2, (count - 1) * pitch + w_bus / 2)
        x_port = x_bus - w_bus / 2 - 1500.0
        feed = box(x_port, y_mid - w_bus / 2, x_bus, y_mid + w_bus / 2)
        shared_port = Point(x_port, y_mid).buffer(r_in, quad_segs=64)
        bus_geo = unary_union(stubs + [bus, feed, shared_port])
        parts.append(bus_geo)
    fluid = unary_union(parts).buffer(0)
    os.makedirs(outdir, exist_ok=True)
    aname = '%s_x%d' % (name, count)
    dxf = os.path.join(outdir, aname + '_UNASSESSED.dxf')   # named by its verdict below (code review 4, F4)
    V = {}
    comps = list(fluid.geoms) if hasattr(fluid, 'geoms') else [fluid]
    V['components'] = dict(n=len(comps), expected=1 if shared else count, ok=len(comps) == (1 if shared else count))
    # In each element, once its outlet stretch is cut, the product port must not share a component with its filtrate
    # port(s) or with the inlet.
    o = arr['outlet']
    leg, gx0 = arr['legs'][-1], arr['gx'][-1]
    xl = o['x_e_local'] - spec.L_out / 2.0
    hw = float(arr['wc'][-1]) / 2.0 + spec.l_g + float(np.nan_to_num(arr['ws'][int(arr['N']) - 2])) + 30.0
    cut0 = L.to_global(box(xl - 5.0, -hw, xl + 5.0, hw), leg, gx0)
    sep = []
    for j in range(count):
        cut = affinity.translate(cut0, 0.0, j * pitch)
        pieces = fluid.difference(cut)
        pcs = list(pieces.geoms) if hasattr(pieces, 'geoms') else [pieces]

        def comp(pt):
            P = Point(pt[0], pt[1] + j * pitch)
            return next((i for i, c in enumerate(pcs) if c.contains(P)), None)
        ic, ia, ib = comp(o['ports']['C']), comp(o['ports']['A']), comp(o['ports']['B'])
        iin = comp((x_in, 0.0))
        ok = None not in (ic, ia, ib, iin) and ic not in (ia, ib, iin) and iin not in (ia, ib)
        sep.append(bool(ok))
    V['product_port_separated_in_every_element'] = dict(ok=all(sep), per_element=sep)
    # The element's 2 mm rule between its inlet port and every outlet port is applied between copies
    # (inter_copy_port_spacing; code review 4 F10; third-party review 2 of rc7, C4).
    o_r = {k: float((o.get('port_r') or {}).get(k, spec.port_r)) for k in o['ports']}
    r_own = float(arr['inlet'].get('port_r', r_in))     # the element's own inlet port as drawn (1020 um, Dinh's filter)
    inl = [((x_port, y_mid), r_in, None)] if shared else [((x_in, j * pitch), r_own, j) for j in range(count)]
    V['inter_copy_port_spacing'] = inter_copy_port_spacing(inl, o['ports'], o_r, count, pitch, L.PORT_GAP_IN_OUT)
    if bus_geo is not None:
        own = unary_union([affinity.translate(port_in, 0.0, j * pitch).buffer(2.0) for j in range(count)])
        stray = bus_geo.intersection(unary_union(els)).difference(own).area
        gap_min = min(bus_geo.difference(own.buffer(w_bus)).distance(e.difference(affinity.translate(port_in, 0.0, j * pitch).buffer(w_bus + 5)))
                      for j, e in enumerate(els))
        V['manifold_touches_elements_only_at_their_ports'] = dict(ok=stray < 1.0 and gap_min >= spec.F_min - 1e-6,
                                                                 stray_area_um2=float(stray), min_clearance_um=float(gap_min))
        # The throughput split takes the element's resistance from its record and the bus and stub resistances from
        # the duct law.
        R_el = rec['design'].get('R_element_SI')
        h = spec.h
        if R_el:
            Rb = S.rr(w_bus, h) * pitch * S.UM
            Rs_ = S.rr(w_bus, h) * (x_in - (x_bus - w_bus / 2)) * S.UM
            n = count
            # The bus has a node at each stub, and the feed enters at the middle node(s).
            A = np.zeros((n, n)); rhs = np.zeros(n)
            g_el = 1.0 / (Rs_ + R_el)
            for i in range(n):
                A[i, i] += g_el
                if i > 0:
                    A[i, i] += 1 / Rb; A[i, i - 1] -= 1 / Rb
                if i < n - 1:
                    A[i, i] += 1 / Rb; A[i, i + 1] -= 1 / Rb
            mid = [n // 2] if n % 2 else [n // 2 - 1, n // 2]
            for m_ in mid:
                rhs[m_] = 1.0 / len(mid)
            P = np.linalg.solve(A, rhs)
            q = P * g_el
            V['throughput_split'] = dict(min_over_mean=float(q.min() / q.mean()), max_over_mean=float(q.max() / q.mean()))
            # The Reynolds gate is evaluated at the busiest element's share of the stated total flow, because the
            # element's Re scales with its flow.
            g_el = rec['validation'].get('gates', {})
            Re_el = g_el.get('Re_max')
            share = count * float(q.max() / q.sum())          # the busiest element's flow over its stated flow
            if Re_el is not None:
                Re_arr = float(Re_el) * share
                V['reynolds_per_element'] = dict(note='total flow = %d x the element\'s stated flow' % count,
                                                 Re_max=Re_arr, ok=bool(Re_arr <= BU.RE_MAX))
            if g_el.get('Re_gap_max') is not None:              # the Stokes gates scale with the flow too (CIF_v2 M14)
                rg, rd = float(g_el['Re_gap_max']) * share, float(g_el['Re_delta_max']) * share
                V['stokes_per_element'] = dict(note='Re_gap and Re_delta at the busiest element\'s flow', Re_gap_max=rg,
                                               Re_delta_max=rd, ok=bool(rg <= BU.RE_STOKES and rd <= BU.RE_STOKES))
    # Copies or a manifold that merge, or a product port joined to a filtrate port, are refused before the array's
    # photomask is written, so the composition is NOT BUILT, as for a pitch below the floor (the clearance guarantee).
    coll_ = [k for k in ('components', 'product_port_separated_in_every_element', 'inter_copy_port_spacing',
                         'manifold_touches_elements_only_at_their_ports') if k in V and not V[k]['ok']]
    if coll_:
        raise ValueError('the composed array\'s parts merge or disconnect, or a copy\'s inlet crowds another copy\'s '
                         'outlet ports (%s): %s; no photomask written'
                         % (', '.join(coll_), '; '.join('%s %s' % (k, {kk: vv for kk, vv in V[k].items() if kk != 'ok'})
                                                       for k in coll_)))
    io_dxf.write_dxf(fluid, dxf)
    if chip:
        cw, ch = chip
        bb = fluid.bounds
        # The chip area is a packaging budget, which is reported and is never a pass condition (App. B thesis; C8 CQ5
        # "reports the result without failing").
        V['fits_chip'] = dict(ok=(bb[2] - bb[0]) <= cw - 2 * margin_um and (bb[3] - bb[1]) <= ch - 2 * margin_um,
                              footprint_mm=[(bb[2] - bb[0]) / 1e3, (bb[3] - bb[1]) / 1e3], reported_only=True)
    if do_drc:
        aj = os.path.join(outdir, aname + '_audit.json')
        if os.path.exists(aj):
            os.remove(aj)
        r = subprocess.run([sys.executable, BU.DRC, aj, dxf, '%.2f' % (spec.F_min - 0.1), '--islands'], capture_output=True, text=True)
        try:
            dr = json.load(open(aj)) if r.returncode == 0 else dict(error=r.stderr[-500:])
        except Exception as e:
            dr = dict(error=str(e))
        V['drc'] = dict(ok=('error' not in dr and dr.get('n_lost_open') == 0 and dr.get('n_lost_solid') == 0),
                        **{k: dr.get(k) for k in ('n_components', 'n_lost_open', 'n_lost_solid', 'd_um', 'error')})
    notes.append('each element keeps its own outlet ports and the off-chip outlet tubing it was validated with '
                 '(R_ext centre %.3g, side %.3g Pa s/m^3); different tubing changes its last rows' % (
                     float(getattr(spec, 'R_ext_c', 0.0)), float(getattr(spec, 'R_ext_s', 0.0))))
    ok = all(v.get('ok', True) for v in V.values() if isinstance(v, dict) and not v.get('reported_only'))
    failed = [k for k, v in V.items() if isinstance(v, dict) and not v.get('ok', True) and not v.get('reported_only')]
    over = V.get('fits_chip') if isinstance(V.get('fits_chip'), dict) and not V['fits_chip']['ok'] else None
    D_el, F_el = rec['validation']['realised_cutoff']['D_target_um'], rec['validation']['network']['F_as_drawn']
    from .pipeline import code_stamp
    code = code_stamp()
    el_code = rec.get('code') or {}
    if el_code.get('sha256') != code['sha256']:
        notes.append('the element was built by code %s; this array is composed by code %s' % (el_code.get('sha256'), code['sha256']))
    a_el = rec.get('assessment', {})
    i_m = verdict.find('[1-D model')
    model_txt = (' ' + verdict[i_m:verdict.index(']', i_m) + 1]) if i_m >= 0 and ']' in verdict[i_m:] else ''
    ar_txt = ''
    if a_el.get('over_aspect_ratio'):                    # the copies carry the element's geometry and its aspect ratio
        i_a = verdict.find(' -- ABOVE THE PROCESS ASPECT-RATIO LIMIT')
        ar_txt = verdict[i_a:] if i_a >= 0 else ' -- ABOVE THE PROCESS ASPECT-RATIO LIMIT (the element\'s; reported)'
        ar_txt = ar_txt.split(' -- OVER THE CHIP AREA')[0]
    if a_el.get('over_chip_area'):
        notes.append('the element is over its own chip area (its record); the array is fitted to its own chip area below')
    out = dict(name=aname, element=name, element_verdict=verdict, element_code=el_code,
               element_physics_audit=(rec.get('validation') or {}).get('physics_audit', {}).get('status'),
               code=code, count=count, pitch_um=pitch,
               dxf=os.path.basename(dxf), dxf_sha256=io_dxf.sha256(dxf),
               inlet='shared' if shared else 'separate',
               bus_width_um=w_bus if shared else None, checks=V, notes=notes,
               over_aspect_ratio=bool(a_el.get('over_aspect_ratio')), over_chip_area=bool(over),
               verdict=('PASS: %d copies of the validated element (D %.4f um, F %.5f on its drawn mask), outlets per element%s%s%s'
                        % (count, D_el, F_el, model_txt, ar_txt,
                           (' -- OVER THE CHIP AREA (%.2f x %.2f mm drawn; reported, not refused or failed)'
                            % tuple(over['footprint_mm'])) if over else '')) if ok else
                       ('FAIL (%s): the array is not a valid composition of the element' % ', '.join(failed)),
               exit_code=0 if ok else 2)
    # The photomask is never under the deliverable name unless the array passed (the clearance guarantee).
    from .pipeline import name_mask
    out['dxf'] = name_mask(outdir, aname, out['exit_code'])
    json.dump(out, open(os.path.join(outdir, aname + '.json'), 'w'), indent=1, default=float)
    # the figure of the whole array
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from opencif.figures import drawing as RD
    bb = fluid.bounds
    fig, ax = plt.subplots(figsize=(12, 12 * (bb[3] - bb[1]) / max(bb[2] - bb[0], 1) + 1))
    RD.draw(ax, fluid, (bb[0] - 300, bb[1] - 300, bb[2] + 300, bb[3] + 300), style='layout')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title('%s: %d x %s, %s inlet; %s' % (aname, count, name, out['inlet'], out['verdict'][:60]), fontsize=9, loc='left')
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, aname + '_layout.png'), dpi=150)
    plt.close(fig)
    return out
