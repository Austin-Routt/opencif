"""test_trace_anchors.py: the micrograph traces (provenance/rc7_traces) re-run from their inputs.

The traces use the rc7 trace code and its window rule, and they cover every shipped mask that the micrographs show,
including the masks of the rc7 presets. The suite has six layers.
  Layer 1  The code. code_stamp() equals the anchors' stamp, because the traces are of the masks that this code draws.
  Layer 2  The inputs. The suite checks the four paper rasters (sha256) and reads every scale bar with its anti-aliased
           ends. It checks Xia's scale (the printed widths of panels c and g, and the curvature of panel d's loop arcs
           against the mask's loop radii) and Fig. 3a's leg spacing against panel d. Gifford's Fig. 6 is measured again
           from its raster (provenance/rc6_scripts/measure_gifford_fig6.py) and compared with the measurement that the
           presets carry, and its turns (provenance/rc7_traces/code/gifford6_turns.py) are compared with the
           measurement that the leg axes of the Fig. 6 presets come from.
  Layer 3  The shipped masks on the micrographs. Each panel is registered again from the raster and the DXF (by
           features and, independently, by chamfer). Every panel must reproduce its anchored registration and metric,
           because the code is deterministic. A claim must also pass its limits, which are V1 (the two registrations
           agree), V2 (the nulls, either the metric's or the chamfer's, as the anchors say), the median and
           90th-percentile limits, and V5 (the scale that the image prefers, along and across the flow). A recorded
           panel is reproduced and its recorded ranges are checked. Xia's Fig. 3a whole-device registration (V11 group
           shifts, V12 the panels' boxes) is re-run and compared.
  Layer 4  The presets' layouts. The Gifford presets carry the layouts measured on their figures
           (provenance/rc6_scripts/gifford_layout_measurement.json). Each variant preset differs from its shipped
           preset in exactly the declared fields and plans to the rows that its record holds.
  Layer 5  The variant masks on the micrographs, as in layer 3. The layer has no panels, because rc6's two Xia
           variants became the shipped xia2016_rs7_aligned and rc7's rejected Fig. 4 reading is figure-only. Gifford
           Fig. 6's one-map trace of the four legs and their turns belongs to layer 3, because the Fig. 6 presets take
           the micrograph's leg axes.
  Layer 6  V13, the taper-direction sign test on every Gifford Fig. 4 reading
           (provenance/rc7_traces/code/taper_sign.py), re-run from the maps of the trace records
           (provenance/rc7_traces/results/gifford/) and the DXFs (the shipped runs/ and the variant builds). Each
           reading must reproduce its anchored verdict, positions, and slopes. The neighbor of the nested reproduction
           must PASS, that of the path-order reading B must FAIL, and leg 4 (the control) must PASS on every mask
           (V13_REQUIRED, held here as well as in the anchors). The layer needs gif5-001.png and the variant builds
           (pathorder, start51, fold20).
Every panel's coverage reproduces its anchor, a claim's coverage must be at least its coverage_min, and its metric
nulls are judged on the censored ratios (censored_null_ratio_min) where its limits say so (A7).

Two inputs are not in the repository, because the papers' figures are the publishers':
  OPENCIF_PAPER_RASTERS   a folder with dinh5-000.png, gif5-001.png, gif7-000.png, and xia5-000.png, made by
                          `pdfimages -png` from the three PDFs (provenance/rc6_traces/README.md gives the commands and
                          the hashes, and rc7_traces uses the same)
  OPENCIF_TRACE_VARIANTS  a folder with the variant builds, <name>/<name>.dxf and .json for each variant (by default
                          provenance/rc7_traces/variants/runs). Layer 5 runs nothing, because no variant panel is
                          anchored, whereas layer 6 reads the DXFs of the three Fig. 4 variants. A variant is built with
    python3 provenance/rc7_traces/code/build_variant.py . provenance/rc7_traces/variants/<name>.yaml <folder>/<name>
A missing input is named and its checks are not run (exit 4 if nothing failed), and a panel missing from the anchors is
a named missing input too.

Usage: python3 tests/test_trace_anchors.py [--layers 1,2,3,4,5,6] [--quick]
The option --quick skips the V5 scale scans and the Fig. 6 re-measurement. The exit codes are
  0   every check passed
  1   a check failed
  4   an input is missing (named), and nothing failed"""
import json
import os
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TR = os.path.join(ROOT, 'provenance', 'rc7_traces')
CODE_DIR = os.path.join(TR, 'code')
# The figure measurements that the traces and presets rest on are rc6's (rc7 re-measured none of them).
MEAS = os.path.join(ROOT, 'provenance', 'rc6_traces', 'measurements')
sys.path[:0] = [CODE_DIR, ROOT]

RESULTS, MISSING = [], []


def check(name, ok, detail=''):
    RESULTS.append((name, bool(ok), detail))
    print('%-6s %s  %s' % ('ok' if ok else 'FAIL', name, detail), flush=True)


def _conv(o):
    """Convert a numpy array, scalar, or bool to its Python value for json.dump, and raise TypeError for anything
    else."""
    import numpy as np
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    raise TypeError(type(o))


def missing(what):
    """Record and print an input that is missing, whose checks are then not run."""
    MISSING.append(what)
    print('MISSING %s' % what, flush=True)


def main():
    """Run the layers that --layers names (all six by default) and exit with 0, 1, or 4."""
    args = sys.argv[1:]
    layers = {1, 2, 3, 4, 5, 6}
    if '--layers' in args:
        layers = {int(v) for v in args[args.index('--layers') + 1].split(',')}
    quick = '--quick' in args
    A = json.load(open(os.path.join(TR, 'trace_anchors.json')))
    import panel_core
    panel_core.DRAW = False
    if 1 in layers:
        from opencif.pipeline import code_stamp
        cs = code_stamp()
        # The anchors were traced on RC8's photomasks. A later release whose evidence shows every photomask byte for
        # byte RC8's (validation/RELEASE.json) traces the same masks, so its own stamp stands for the anchors' too.
        rel_p = os.path.join(ROOT, 'validation', 'RELEASE.json')
        rel = json.load(open(rel_p)) if os.path.exists(rel_p) else {}
        cnt = rel.get('counts') or {}
        same_masks = bool(rel.get('code_stamp') == cs['sha256'] and rel.get('rc8_code_stamp') == A['code_stamp']['sha256']
                          and cnt.get('photomasks') and cnt.get('identical_to_rc8') == cnt.get('photomasks'))
        check('L1 code stamp', cs['sha256'] == A['code_stamp']['sha256'] or same_masks,
              '%s (anchors %s%s)' % (cs['sha256'], A['code_stamp']['sha256'],
                                     '; the release draws RC8\'s photomasks byte for byte' if same_masks else ''))
    ras = os.environ.get('OPENCIF_PAPER_RASTERS', '')
    have_ras = bool(ras) and all(os.path.exists(os.path.join(ras, f)) for f in A['rasters'])
    if (layers & {2, 3, 5, 6}) and not have_ras:         # layer 6 reads gif5-001.png (RC8 code review, m1)
        missing('OPENCIF_PAPER_RASTERS with %s' % ', '.join(A['rasters']))
    if 2 in layers and have_ras:
        run_layer2(A, ras, quick)
    if 3 in layers and have_ras:
        run_traces(A, ras, quick, 'shipped')
    if 4 in layers:
        run_layer4(A)
    if 5 in layers and have_ras:
        run_traces(A, ras, quick, 'variants')
    if 6 in layers and have_ras:
        run_layer6(A, ras)
    n_fail = sum(1 for r in RESULTS if not r[1])
    print('\n%d checks, %d failed, %d inputs missing' % (len(RESULTS), n_fail, len(MISSING)))
    sys.exit(1 if n_fail else (4 if MISSING else 0))


# ------------------------------------------------------------------------------------------------------------ layer 2
def run_layer2(A, ras, quick):
    """Check the inputs of layer 2, which are the rasters, the scale bars, Xia's scale and leg spacing, and Gifford's
    Fig. 6 and its turns measured again (not with quick)."""
    import trace_lib as T
    for f, d in A['rasters'].items():
        h = T.sha256(os.path.join(ras, f))
        check('L2 raster %s' % f, h == d['sha256'], h[:16])
    im = {f: T.Raster(os.path.join(ras, f)).im for f in A['rasters']}
    B = A['bars']
    for rr in B['dinh_fig2B']['rows']:
        b = T.bar_length_aa(im['dinh5-000.png'], rr[0], rr[1], *B['dinh_fig2B']['cols'])
        check('L2 Dinh Fig. 2B bar rows %d-%d' % tuple(rr), b is not None and abs(b[0] - B['dinh_fig2B']['px_aa']) <= 0.05
              and b[3] == B['dinh_fig2B']['full_px'], '%.3f px (full pixels %d)' % (b[0], b[3]) if b else 'not found')
    g = B['gifford_fig4_inset']
    b = T.bar_length_aa(im['gif5-001.png'], g['rows'][0], g['rows'][1], *g['cols'])
    check('L2 Gifford Fig. 4 inset bar', b is not None and abs(b[0] - g['px_aa']) <= 0.05, '%.3f px' % b[0] if b else 'not found')
    x = B['xia_fig3_panel_i']
    b = T.bar_length_aa(im['xia5-000.png'], x['rows'][0], x['rows'][1], *x['cols'], dark=False)
    check('L2 Xia panel-i bar (not the scale of panels b-g)', b is not None and abs(b[0] - x['px_aa']) <= 0.05,
          '%.3f px (%.4f µm/px)' % (b[0], 250.0 / b[0]) if b else 'not found')
    x3 = B['xia_fig3a']
    b = T.bar_length_aa(im['xia5-000.png'], x3['rows'][0], x3['rows'][1], *x3['cols'])
    check('L2 Xia Fig. 3a bar', b is not None and abs(b[0] - x3['px_aa']) <= 0.05, '%.3f px' % b[0] if b else 'not found')
    S = A['scales']
    import panels_xia as PX
    ev = PX.scale_evidence(im['xia5-000.png'])
    check('L2 Xia printed-width scale (panels c and g)', abs(ev['um_per_px_printed'] - S['xia']['printed_widths_um_per_px']) <= 0.002,
          '%.4f µm/px (anchored %.4f)' % (ev['um_per_px_printed'], S['xia']['printed_widths_um_per_px']))
    sep = ev['panel_d_axis_separation_px']
    check('L2 Xia leg spacing (Fig. 3a 901 ± 28 µm) at the rc6 scale', sep is not None and abs(sep * S['xia']['rc6_um_per_px'] - 901.0) <= 28.0,
          '%.1f µm' % (sep * S['xia']['rc6_um_per_px']) if sep else 'not found')
    xd = os.path.join(ROOT, 'runs', 'xia2016_rs7_depicted', 'xia2016_rs7_depicted.dxf')
    if os.path.exists(xd):
        stem = os.path.join(tempfile.mkdtemp(prefix='arcs_'), 'arcs')
        subprocess.run([sys.executable, '-I', os.path.join(CODE_DIR, 'xia_arc_scale.py'), os.path.join(ras, 'xia5-000.png'), xd, stem],
                       check=True, capture_output=True)
        a = json.load(open(stem + '.json'))['scale_um_per_px']
        check('L2 Xia scale from the loop arcs\' curvature (panel d vs the mask\'s loop radii)',
              abs(a['mean'] - S['xia']['arc_um_per_px']) <= 0.005 and a['sd'] <= 0.01,
              '%.4f ± %.4f µm/px over %d arcs (anchored %.4f)' % (a['mean'], a['sd'], a['n'], S['xia']['arc_um_per_px']))
    else:
        missing('runs/xia2016_rs7_depicted/xia2016_rs7_depicted.dxf')
    check('L2 Xia rc6 scale = the mean of its three determinations (data/r6/rc6_values.json)',
          abs(S['xia']['rc6_um_per_px'] - json.load(open(os.path.join(ROOT, 'opencif', 'data', 'r6', 'rc6_values.json')))['xia']['scale_um_per_px']) < 1e-9,
          '%.4f µm/px' % S['xia']['rc6_um_per_px'])
    if quick:
        print('   (--quick: Gifford Fig. 6 not measured again)')
        return
    out = os.path.join(tempfile.mkdtemp(prefix='g6m_'), 'm.json')
    t0 = time.time()
    subprocess.run([sys.executable, '-I', os.path.join(ROOT, 'provenance', 'rc6_scripts', 'measure_gifford_fig6.py'),
                    os.path.join(ras, 'gif7-000.png'), out], check=True, capture_output=True)
    m, m0 = json.load(open(out)), json.load(open(os.path.join(ROOT, 'provenance', 'rc6_scripts', 'gifford_fig6_measurement.json')))
    same = all(abs(m[k] - m0[k]) <= 1e-6 * max(1.0, abs(m0[k])) for k in ('um_per_px', 'row_pitch_um', 'l_g_um', 'rows_direction_deg'))
    check('L2 Gifford Fig. 6 measured again from its raster = the presets\' measurement (%.0f s)' % (time.time() - t0), same,
          'pitch %.3f µm, l_g %.2f µm, %.4f µm/px, rows %.3f°' % (m['row_pitch_um'], m['l_g_um'], m['um_per_px'], m['rows_direction_deg']))
    check('L2 Gifford Fig. 6: the 250 µm bar (cut by the raster) bounds the scale from above',
          m['scale_bar']['runs_off_the_right_edge'] and m['um_per_px'] <= m['scale_bar']['um_per_px_upper_bound'],
          '%.4f <= %.4f µm/px' % (m['um_per_px'], m['scale_bar']['um_per_px_upper_bound']))
    # The turn walls are followed again and fitted about one center, and the fit is compared with
    # measurements/gifford_fig6_turns.json, from which the leg axes of the Fig. 6 presets (legs 2, 4, 6) come (after the
    # 2026-10-17 review).
    out = os.path.join(tempfile.mkdtemp(prefix='g6t_'), 't.json')
    subprocess.run([sys.executable, '-I', os.path.join(CODE_DIR, 'gifford6_turns.py'), ROOT, ras, out], check=True, capture_output=True)
    t, t0_ = json.load(open(out)), json.load(open(os.path.join(MEAS, 'gifford_fig6_turns.json')))
    same = (t['walls_in_the_turns'] == t0_['walls_in_the_turns'] and
            max(abs(a - b) for a, b in zip(t['centre_uv_px'], t0_['centre_uv_px'])) <= 1e-6 and t['leg_Y_um'] == t0_['leg_Y_um'])
    check('L2 Gifford Fig. 6 turns measured again: %d walls about one centre (residual SD %.2f px), %.1f µm beyond leg 1; '
          'the presets\' leg axes follow' % (t['walls_in_the_turns'], t['residual_sd_px'], t['centre_from_leg1_axis_um']),
          same and t['residual_sd_px'] <= 0.25, 'leg_Y %s' % t['leg_Y_um'])


# ------------------------------------------------------------------------------------------------------- layers 3, 5
def variant_paths(name):
    """Return the preset, DXF, and record paths of the variant build name, or None when the build is missing."""
    base = os.environ.get('OPENCIF_TRACE_VARIANTS', os.path.join(TR, 'variants', 'runs'))
    d = os.path.join(base, name)
    p = dict(preset=os.path.join(TR, 'variants', name + '.yaml'), dxf=os.path.join(d, name + '.dxf'), record=os.path.join(d, name + '.json'))
    return p if os.path.exists(p['dxf']) and os.path.exists(p['record']) else None


def judge(rec, A):
    """Check that the trace record rec reproduces its panel's anchored values, and then its claim limits or its
    recorded ranges."""
    import numpy as np
    name = rec['name']
    P = A['panels'].get(name)
    if P is None:
        missing('trace_anchors.json has no entry for the panel %s' % name)
        return
    a, lim = P['anchored'], P['limits']
    v1 = rec['V1_feature_vs_chamfer_px']
    rep = (abs(v1['along'] - a['V1_along_px']) <= 0.05 and abs(v1['across'] - a['V1_across_px']) <= 0.05
           and abs(rec['metric']['median_abs_um'] - a['median_abs_um']) <= 0.03
           and abs(rec['V2_chamfer']['best_px'] - a['chamfer_best_px']) <= 0.01)
    check('L%s %s reproduces its anchors' % (P['layer'], name), rep,
          'V1 %.2f / %.2f px, median %.3f µm, chamfer %.3f px (anchored %.2f / %.2f, %.3f, %.3f)' % (
              v1['along'], v1['across'], rec['metric']['median_abs_um'], rec['V2_chamfer']['best_px'],
              a['V1_along_px'], a['V1_across_px'], a['median_abs_um'], a['chamfer_best_px']))
    tag = 'claim' if P['role'] == 'claim' else 'recorded'
    wr = rec.get('window_rule')
    if isinstance(wr, dict):                         # the metric half-window at most half the drawn gap
        check('L%s %s window rule' % (P['layer'], name), wr['ok'],
              'half-window %.2f µm, half the gap %.2f µm' % (rec['metric_window_um'], wr['half_gap_um']))
    if 'v1' in lim:
        check('L%s %s (%s) V1 the two registrations %s' % (P['layer'], name, tag, 'agree' if lim['v1'] else 'disagree, as recorded'),
              v1['agree'] == lim['v1'], 'along %.2f, across %.2f px' % (v1['along'], v1['across']))
    med, p90 = rec['metric']['median_abs_um'], rec['metric']['p90_abs_um']
    if 'median_max_um' in lim:
        check('L%s %s (%s) median / p90' % (P['layer'], name, tag), med <= lim['median_max_um'] and p90 <= lim['p90_max_um'],
              '%.2f / %.2f µm (limits %.1f / %.1f)' % (med, p90, lim['median_max_um'], lim['p90_max_um']))
    if 'metric_null_ratio_min' in lim:
        nh = rec['metric_null']['half_pitch_along']['median_abs_um']
        na = rec['metric_null']['5um_across']['median_abs_um']
        r = lim['metric_null_ratio_min']
        check('L%s %s (%s) V2 metric nulls' % (P['layer'], name, tag), min(nh, na) >= r * med,
              'median %.2f µm; half a pitch %.2f, 5 µm across %.2f (both >= %.2f x)' % (med, nh, na, r))
    if 'censored_null_ratio_min' in lim:              # the nulls over every eligible sample (A7)
        cr = rec['metric_null'].get('censored_ratios') or {}
        r = lim['censored_null_ratio_min']
        ok_ = (cr.get('half_pitch_along') is not None and cr.get('5um_across') is not None   # either may be None (n5)
               and min(cr['half_pitch_along'], cr['5um_across']) >= r)
        check('L%s %s (%s) V2 metric nulls, censored' % (P['layer'], name, tag), ok_,
              'half a pitch %s, 5 µm across %s (both >= %.2f x)' % (cr.get('half_pitch_along'), cr.get('5um_across'), r))
    cov = rec['metric'].get('coverage')
    if a.get('coverage') is not None:                 # the coverage anchor (A7)
        check('L%s %s coverage reproduces its anchor' % (P['layer'], name), cov is not None and abs(cov - a['coverage']) <= 0.002,
              '%s (anchored %.4f)' % (('%.4f' % cov) if cov is not None else 'none', a['coverage']))
    if 'coverage_min' in lim:
        check('L%s %s (%s) coverage' % (P['layer'], name, tag), cov is not None and cov >= lim['coverage_min'],
              '%s of the eligible samples found an outline (>= %.2f)' % (('%.3f' % cov) if cov is not None else 'none', lim['coverage_min']))
    if 'chamfer_null_ratio_min' in lim:
        c = rec['V2_chamfer']
        worst = min(min(c['half_pitch_along_px']), min(c['five_um_across_px']))
        r = lim['chamfer_null_ratio_min']
        check('L%s %s (%s) V2 chamfer nulls' % (P['layer'], name, tag), worst >= r * c['best_px'],
              'best %.3f px; half a pitch %s, 5 µm across %s (both >= %.2f x)' % (c['best_px'], np.round(c['half_pitch_along_px'], 3),
                                                                                  np.round(c['five_um_across_px'], 3), r))
    if 'V5_scale_scan' in rec and 'v5_along' in lim:
        v = rec['V5_scale_scan']
        (la, ha), (lc, hc) = lim['v5_along'], lim['v5_across']
        check('L%s %s (%s) V5 scale%s' % (P['layer'], name, tag, (' (%s)' % lim['v5_note']) if lim.get('v5_note') else ''),
              la <= v['best_along'] <= ha and lc <= v['best_across'] <= hc,
              'along %.4f (%.3f-%.3f), across %.4f (%.3f-%.3f)' % (v['best_along'], la, ha, v['best_across'], lc, hc))
    for key, (lo, hi, path, what) in (lim.get('measured') or {}).items():
        val = rec
        for p_ in path:
            val = val[p_] if not isinstance(val, list) else val[int(p_)]
        check('L%s %s (%s) %s' % (P['layer'], name, tag, what), lo <= val <= hi, '%.3f (%.3f-%.3f)' % (val, lo, hi))


def run_traces(A, ras, quick, which):
    """Re-run the traces of layer 3 (which is 'shipped') or layer 5 (which is 'variants') and judge every panel
    record."""
    import trace_lib as T
    scan = not quick
    tmp = tempfile.mkdtemp(prefix='traces_')
    os.environ['XIA_RASTERS'] = ras
    from opencif import spec as SP
    L = '3' if which == 'shipped' else '5'
    jobs = A['runs'][which]
    for job in jobs:
        t0 = time.time()
        kind = job['kind']
        paths = None
        if which == 'shipped':
            n = job['mask']
            paths = dict(preset=os.path.join(ROOT, 'presets', n + '.yaml'), dxf=os.path.join(ROOT, 'runs', n, n + '.dxf'),
                         record=os.path.join(ROOT, 'runs', n, n + '.json'))
            if not os.path.exists(paths['dxf']):
                missing('runs/%s/%s.dxf (build the preset)' % (n, n))
                continue
        else:
            paths = variant_paths(job['mask'])
            if paths is None:
                missing('the variant build %s (OPENCIF_TRACE_VARIANTS)' % job['mask'])
                continue
        label = job['label']
        tag = job['tag']
        recs = []
        kw_ = dict(paths, label=label)                   # every run names its mask (the shipped runs included)
        if kind == 'dinh':
            import panels_dinh as PD
            import panels_dinh_v2 as PV
            ctx = PD.setup(ROOT, ras, tag=tag, **kw_)
            r4 = PV.panel4(ctx, tmp, scan)
            recs = [r4, PV.panel3(ctx, tmp, r4, scan), PV.panel2(ctx, tmp, scan), PV.panel1(ctx, tmp, scan)]
        elif kind == 'gifford4':
            import panels_gifford as PG
            recs = [PG.run(ROOT, ras, tmp, tag=tag, scan=scan, **kw_)]
        elif kind == 'gifford6':
            import panels_gifford6 as P6
            recs = P6.run(ROOT, ras, tmp, tag=tag, scan=scan, one_map=job.get('one_map', False), turns=job.get('turns', False), draw=False,
                          **kw_)
            m = recs[0]['measured']
            ec = m['end_check']
            check('L%s Gifford Fig. 6 (%s): END_CHECK -- on every row the last pill (image) and the last post (mask) stand one gap '
                  'before the divider wall, no further post beyond (the along-flow pairing)' % (L, tag), ec['ok'],
                  'image gaps %s µm, mask %s µm' % (' / '.join('%.1f' % r['image_gap_um'] for r in ec['rows'] if r['image_gap_um']),
                                                   sorted({round(r['mask_gap_um'], 2) for r in ec['rows'] if r['mask_gap_um']})))
            check('L%s Gifford Fig. 6 (%s): the legs identified as 7, 5, 3, 1 (flow direction and side widths)' % (L, tag),
                  m['identification']['best']['legs'] == [7, 5, 3, 1] and not m['identification']['best']['mirrored'],
                  'offset SD %.2f µm (mirrored order: %.1f)' % (m['identification']['best']['offset_sd_um'],
                                                                 [c['offset_sd_um'] for c in m['identification']['candidates']
                                                                  if c['mirrored'] and c['same_flow_direction']][0]))
        elif kind == 'xia':
            import panels_xia as PX
            import panels_xia_more as PM
            rec0 = json.load(open(paths['record']))
            mask = T.Mask(paths['dxf'])
            if mask.sha256 != rec0['files']['dxf_sha256']:
                check('L%s %s: the DXF is its record\'s' % (L, job['mask']), False, mask.sha256[:16])
                continue
            sp, op = SP.to_engine(SP.load(paths['preset']))
            s6 = json.load(open(os.path.join(ROOT, 'opencif', 'data', 'r6', 'rc6_values.json')))['xia']['scale_um_per_px']
            note = ('the rc6 scale %.3f µm/px: the mean of the printed widths, the post columns\' separations and the loop '
                    'curvature (data/r6/rc6_values.json)' % s6)
            f3 = os.path.join(tmp, 'f3a_%s.json' % tag)
            subprocess.run([sys.executable, '-I', os.path.join(CODE_DIR, 'xia_fig3a_device.py'), os.path.join(ras, 'xia5-000.png'),
                            paths['dxf'], paths['record'], f3, f3[:-5] + '.png'], check=True, capture_output=True)
            F3 = json.load(open(f3))
            g = PX.panel_g(mask, sp, rec0, tmp, 1.0 / s6, tag, label, note, scan=scan)
            recs.append(g)
            json.dump(g, open(os.path.join(tmp, 'g_%s.json' % tag), 'w'), default=_conv)
            traced = ['%s:-:g' % os.path.join(tmp, 'g_%s.json' % tag)]
            for letter in 'cdef':
                r_ = PM.panel(letter, mask, sp, rec0, tmp, 1.0 / s6, tag, label, note, ras, scan=scan, fig3a=F3)
                recs.append(r_)
                pth = os.path.join(tmp, '%s_%s.json' % (letter, tag))
                json.dump(r_, open(pth, 'w'), default=_conv)
                traced.append('%s:-:%s' % (pth, letter))
            subprocess.run([sys.executable, '-I', os.path.join(CODE_DIR, 'xia_fig3a_device.py'), os.path.join(ras, 'xia5-000.png'),
                            paths['dxf'], paths['record'], f3, f3[:-5] + '.png'] + traced, check=True, capture_output=True)
            judge_fig3a(json.load(open(f3)), A['fig3a'][tag], L, tag)
        for r in recs:
            judge(r, A)
        print('   (%s %s: %.0f s)' % (kind, tag, time.time() - t0), flush=True)


def judge_fig3a(F, a, L, tag):
    """Check that the whole-device pose of Xia's Fig. 3a and each group's shift (V11) reproduce the anchors, and that
    the boxes (V12) are as anchored."""
    ok = abs(F['pose_whole']['rotation_deg'] - a['pose_whole_rotation_deg']) <= 0.002 and abs(F['pose_whole']['chamfer_px'] - a['pose_whole_chamfer_px']) <= 0.01
    check('L%s Xia Fig. 3a (%s) whole-device pose reproduces its anchors' % (L, tag), ok,
          'rotation %.4f°, chamfer %.3f px (anchored %.4f, %.3f)' % (F['pose_whole']['rotation_deg'], F['pose_whole']['chamfer_px'],
                                                                      a['pose_whole_rotation_deg'], a['pose_whole_chamfer_px']))
    for g, s in a['group_shifts'].items():
        f = F['V11_group_shifts'][g]
        ok = abs(f['shift_across_um'] - s['across_um']) <= 1.0 and (s['along_um'] is None or abs(f['shift_along_um'] - s['along_um']) <= 2.0)
        check('L%s Xia Fig. 3a (%s) V11 %s' % (L, tag, g), ok, 'along %s, across %+.1f µm (anchored %s, %+.1f)' % (
            'invariant' if f['shift_along_um'] is None else '%+.1f' % f['shift_along_um'], f['shift_across_um'],
            'invariant' if s['along_um'] is None else '%+.1f' % s['along_um'], s['across_um']))
    for letter, want in a['V12_ok'].items():
        v = F['V12_boxes'].get(letter, {})
        check('L%s Xia Fig. 3a (%s) V12 panel %s inside its box: %s as anchored' % (L, tag, letter, want), v.get('ok') == want,
              '%.2f of the half-size' % v.get('offset_over_half_size', float('nan')))


# ------------------------------------------------------------------------------------------------------------ layer 6
V13_REQUIRED = dict([('aligned_neighbour', 'PASS'), ('pathorder_neighbour', 'FAIL')] +
                    [('%s_leg4' % t_, 'PASS') for t_ in ('rc7', 'aligned', 'pathorder', 'start51', 'fold20')])


def run_layer6(A, ras):
    """Re-run V13 from each Fig. 4 trace record's map and its DXF, and check the verdict, positions, and slopes
    against the anchors. The decisive verdicts are this test's own (V13_REQUIRED) as well as the anchors' (RC8 code
    review, m2), so the neighbor of the nested reproduction PASSES, that of the path-order reading B FAILS, and leg 4
    (the control) PASSES on every mask."""
    V = A.get('v13') or {}
    if not V:
        missing('trace_anchors.json has no v13 entries (RC8)')
        return
    for key, want in sorted(V13_REQUIRED.items()):
        check('L6 V13 %s is anchored with its required verdict %s' % (key, want),
              key in V and V[key].get('required') == want and V[key].get('verdict') == want,
              ('%s, required %s' % (V[key].get('verdict'), V[key].get('required'))) if key in V else 'absent')
    tmp = tempfile.mkdtemp(prefix='v13_')
    for key, v in sorted(V.items()):
        trace = os.path.join(TR, 'results', v['trace'])
        if v['where'] == 'shipped':
            dxf = os.path.join(ROOT, 'runs', v['mask'], v['mask'] + '.dxf')
        else:
            p_ = variant_paths(v['mask'])
            dxf = p_['dxf'] if p_ else os.path.join(TR, 'variants', 'runs', v['mask'], v['mask'] + '.dxf')
        if not (os.path.exists(trace) and os.path.exists(dxf)):
            missing('V13 %s: %s and %s' % (key, trace, dxf))
            continue
        out = os.path.join(tmp, key + '.json')
        cmd = [sys.executable, '-I', os.path.join(CODE_DIR, 'taper_sign.py'), os.path.join(ras, 'gif5-001.png'), trace, dxf,
               v['leg'], key, out, '%g' % v['window_um'], '5', '3'] + ([v['x_min']] if v.get('x_min') else [])
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode or not os.path.exists(out):
            check('L6 V13 %s runs' % key, False, res.stderr[-200:])
            continue
        r = json.load(open(out))
        rep = (r['verdict'] == v['verdict'] and r['positions'] == v['positions'] and
               (v['slope_image'] is None or abs(r['slope_image'] - v['slope_image']) <= 1e-9) and
               (v['slope_mask'] is None or abs(r['slope_mask'] - v['slope_mask']) <= 1e-9))
        check('L6 V13 %s reproduces its anchor (%s)' % (key, v['verdict']), rep,
              '%s, %d positions, image slope %s, mask slope %s' % (r['verdict'], r['positions'], r.get('slope_image'), r.get('slope_mask')))
        want = V13_REQUIRED.get(key) or v.get('required')
        if want:
            check('L6 V13 %s: %s, as required' % (key, want), r['verdict'] == want, r['verdict'])


# ------------------------------------------------------------------------------------------------------------ layer 4
def run_layer4(A):
    """Check the presets' layouts against the figure measurements, and the variant presets against their shipped
    presets (layer 4)."""
    import yaml
    from opencif import spec as SP
    from opencif.physics import schedule as CS
    # The Gifford presets' layouts are the ones their figures show, so the measurement
    # (provenance/rc6_scripts/gifford_layout_measurement.json, from the traces' image-side readings) and the presets
    # agree (after the 2026-10-17 review).
    M = json.load(open(os.path.join(ROOT, 'provenance', 'rc6_scripts', 'gifford_layout_measurement.json')))
    t = json.load(open(os.path.join(MEAS, 'gifford_fig6_turns.json')))
    for n in ('gifford2014_fig6_design', 'gifford2014_fig6_as_fabricated_15um', 'gifford2014_fig6_aligned'):
        r_ = (yaml.safe_load(open(os.path.join(ROOT, 'presets', n + '.yaml'))).get('reproduction') or {})
        ok = ([round(float(v), 1) for v in r_.get('leg_Y', [])] == [round(float(v), 1) for v in M['fig6']['leg_Y_um']] == t['leg_Y_um']
              and r_.get('pc_last_short') == round(M['fig6']['last_leg_end_short_um']))
        check('L4 %s takes the micrograph\'s leg axes and the last leg\'s early end' % n, ok,
              'leg_Y %s, pc_last_short %s' % (r_.get('leg_Y'), r_.get('pc_last_short')))
    for n in ('gifford2014_fig4_design', 'gifford2014_fig4_as_fabricated_20um', 'gifford2014_fig4_as_fabricated_8um',
              'gifford2014_fig4_aligned'):
        lay = yaml.safe_load(open(os.path.join(ROOT, 'presets', n + '.yaml'))).get('layout') or {}
        check('L4 %s takes the inset\'s leg spacing' % n, lay.get('leg_spacing_um') == round(M['fig4']['leg_spacing_um']),
              'leg_spacing_um %s (inset %.1f)' % (lay.get('leg_spacing_um'), M['fig4']['leg_spacing_um']))
    for name, v in A['variants'].items():
        pv = os.path.join(TR, 'variants', name + '.yaml')
        pb = os.path.join(ROOT, 'presets', v['base'] + '.yaml')
        a, b = yaml.safe_load(open(pv)), yaml.safe_load(open(pb))
        diffs = {}
        for sec in sorted(set(a) | set(b)):
            if sec == 'name':
                continue
            if isinstance(a.get(sec), dict) and isinstance(b.get(sec), dict):
                for k in sorted(set(a[sec]) | set(b[sec])):
                    if a[sec].get(k) != b[sec].get(k):
                        diffs['%s.%s' % (sec, k)] = a[sec].get(k)
            elif a.get(sec) != b.get(sec):
                diffs[sec] = a.get(sec)
        check('L4 %s differs from %s only in the declared fields' % (name, v['base']), diffs == v['declared'], json.dumps(diffs))
        sc = CS.plan(SP.to_engine(SP.load(pv))[0])
        check('L4 %s plans to its record\'s rows' % name, sc.N == v['N'] and sc.NM == v['NM'],
              'N %d NM %d (record %d / %d)' % (sc.N, sc.NM, v['N'], v['NM']))


if __name__ == '__main__':
    main()
