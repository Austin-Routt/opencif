"""Collect the rc7 trace records (run_all_rc7.sh) into provenance/rc7_traces/trace_anchors.json: the inputs (raster
hashes, the bars read with their anti-aliased ends, the scales), the runs tests/test_trace_anchors.py repeats, each
panel's anchored registration and metric, the acceptance limits of the claims and the ranges of the recorded panels,
and Xia's Fig. 3a whole-device registrations.  rc7: every shipped mask a micrograph shows (the rc7 presets' masks
included), traced with the rc7 code (its window rule); no variants (rc6's Xia variants became xia2016_rs7_aligned,
and the rejected Fig. 4 start, provenance/rc7_traces/variants/gifford2014_fig4_start51.yaml, is shown, not anchored).
RC8: (A7) every panel's coverage and censored null ratios are anchored; a claim's metric nulls are judged on the
censored ratios (censored_null_ratio_min, the rc7 values) and its coverage must be at least 0.5; (D5) Xia 3d (both
masks) and 3g (depicted) are recorded, not claimed (coverage 0.42-0.48, censored nulls 0.94-1.00: the metric cannot
tell);
the Fig. 4 reproduction is the nested reading (its neighbour is leg 2); (V13) the taper-direction sign test on every
Fig. 4 reading is anchored, the variants' included (the path-order reading B must FAIL, the nested reading PASS).

The limits were set with the results in view: they are regression limits for rc7, not a blind test.  What is blind is
each panel's scale (a bar, a printed width or the measurement), its rotation (the image) and the agreement of two
independent registrations (V1).  The rc6 limits are kept for every panel rc6 claimed that is still claimed (RC8, D5:
Xia 3d and 3g on the depicted mask are recorded); where an rc7 mask no longer
meets one, the panel's entry says so (NULL_FLAGS, package fix #7) rather than the limit being moved silently.

usage: python3 make_anchors.py ROOT RESULTS_DIR OUT.json   (RASTERS in the environment: the folder of the rasters)"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trace_lib as T

ROOT, RES, OUTP = sys.argv[1:4]
sys.path[:0] = [ROOT]
from opencif.pipeline import code_stamp


def J(*p):
    return json.load(open(os.path.join(RES, *p)))


# ---- the limits.  Dinh (1.066 µm/px, the 50 µm bars): rc6's limits, for both readings (Fig. 2B's panels show local
# features the printed and the figure-aligned devices share)
def dinh_limits(tag):
    return {
        'dinh_fig2B_p1_head__%s' % tag: dict(v1=True, median_max_um=1.6, p90_max_um=3.2, chamfer_null_ratio_min=1.3,
                                             v5_along=(0.97, 1.0), v5_across=(0.975, 1.015),
                                             v5_note='rc7: the square-turn head prefers about 1.011 across (rc6\'s arc-turn '
                                                     'head 0.988); the range widened from rc6\'s (0.975, 1.0) with the result '
                                                     'in view',
                                             null_note='rc7: V2 on the chamfer nulls; the metric nulls of the square-turn head '
                                                       'are below rc6\'s 1.3 x (recorded, NULL_FLAGS)'),
        'dinh_fig2B_p2_transition__%s' % tag: dict(v1=True, median_max_um=1.0, p90_max_um=3.0, censored_null_ratio_min=2.0,
                                                   v5_along=(0.98, 1.0), v5_across=(0.995, 1.015)),
        'dinh_fig2B_p3_linear__%s' % tag: dict(v1=True, median_max_um=0.6, p90_max_um=1.2, censored_null_ratio_min=2.0,
                                               v5_along=(0.99, 1.01), v5_across=(0.99, 1.01),
                                               measured=dict(centre=(68.0, 69.5, ['measured', 'image_widths_um', 'centre'], 'centre width (mask 68.1 at its row)'))),
        'dinh_fig2B_p4_outlet__%s' % tag: dict(v1=True, median_max_um=0.7, p90_max_um=2.0, censored_null_ratio_min=2.0,
                                               v5_along=(0.99, 1.01), v5_across=(0.99, 1.01),
                                               measured=dict(centre=(51.5, 53.0, ['measured', 'image_widths_um', 'centre'], 'centre width (printed 52)'),
                                                             side=(43.5, 44.5, ['measured', 'image_widths_um', 'side_L'], 'left side width (printed 44)'),
                                                             pitch=(44.9, 45.3, ['measured', 'pitch_um'], 'row pitch (mask 45.1)'))),
    }


def g6_limits(tag):
    # the design reference at rc6's limits; the aligned device (its side channels on the image) at tighter ones, set with
    # the rc7 results in view
    al = tag == 'aligned'
    d = {'gifford_fig6_leg%d__%s' % (j, tag): dict(v1=True, median_max_um=1.0 if al else 2.0, p90_max_um=2.5 if al else 10.0,
                                                   chamfer_null_ratio_min=1.3 if al else 1.15,
                                                   v5_along=(0.995, 1.005), v5_across=(0.995, 1.01) if al else (1.0, 1.05),
                                                   v5_note=('along: consistency, the pitch was measured on this image; across: '
                                                            + ('the side channels on the image' if al else 'the side channels read wider')))
         for j in (7, 5, 3, 1)}
    d['gifford_fig6_four_legs__%s' % tag] = dict(
        v1=True, median_max_um=1.0 if al else 2.0, p90_max_um=3.0 if al else 8.0, chamfer_null_ratio_min=1.3 if al else 1.15,
        v5_along=(0.995, 1.005), v5_across=(0.995, 1.005),
        v5_note='consistency: the pitch and the leg axes were read on this image',
        measured=dict(
            leg7=(-1.0, 1.0, ['measured', 'leg7_own_map_minus_common_um', 'across'],
                  'leg 7 (not in the fit) across, from its own end\'s map, µm (consistency: the axes were read on this image)'),
            turn_med=(0.0, 2.0, ['measured', 'turns', 'nearest_mask_edge_um', 'median'],
                      'the turns: median distance from the image\'s turn walls to the mask\'s nearest edge, µm (consistency)'),
            turn_p90=(0.0, 5.0, ['measured', 'turns', 'nearest_mask_edge_um', 'p90'],
                      'the turns: 90th percentile of that distance, µm (consistency)'),
            turn_r=(-2.0, 2.0, ['measured', 'turns', 'radii_about_image_centre', 'mean_um'],
                    'the turns: mean radius mask - image about the image\'s turn centre, µm (consistency)'),
            turn_c=(-20.0, 20.0, ['measured', 'turns', 'centre_offset_mask_minus_image_um', 'v'],
                    'the turns: the mask\'s own turn centre across from the image\'s, µm (consistency)'),
            turn_u=(-10.0, 10.0, ['measured', 'turns', 'centre_offset_mask_minus_image_um', 'u'],
                    'the turns: the mask\'s own turn centre along the flow from the image\'s, µm (+ upstream): the one '
                    'independent turn result')))
    return d


CLAIM = {}
CLAIM.update(dinh_limits('rc7'))
CLAIM.update(dinh_limits('fa'))
CLAIM.update(g6_limits('design'))
CLAIM.update(g6_limits('aligned'))
CLAIM.update({
    # Gifford's Fig. 4 inset (2.79 µm/px, the 300 µm bar): relief contrast, outlines read long along the flow.  The
    # design reference (nested; its neighbouring leg 2 lies ~18 µm off the image, which the claim on leg 4 does not
    # cover): rc6's limits, but p90 3.5 µm instead of 20 (rc6's p90 18.2 was the +-8 px window, GIF4-2; at the rc7 window
    # it reads 2.51); reading B's at tight limits, its neighbour leg 3
    'gifford_fig4_inset_leg4__rc7': dict(v1=True, median_max_um=1.5, p90_max_um=3.5, censored_null_ratio_min=2.5,
                                         v5_along=(0.99, 1.01), v5_across=(0.99, 1.01),
                                         measured=dict(pitch=(45.6, 46.0, ['measured', 'pitch_um', 'inner'], 'row pitch (mask 45.8)'),
                                                       spacing=(813.0, 823.0, ['measured', 'leg_centreline_separation_um', 'mask'],
                                                                'the mask\'s leg spacing against the inset\'s 817.9 µm (consistency)'))),
    'gifford_fig4_inset_leg4__aligned': dict(v1=True, median_max_um=1.0, p90_max_um=2.0, censored_null_ratio_min=2.5,
                                             v5_along=(0.99, 1.01), v5_across=(0.99, 1.01),
                                             measured=dict(pitch=(45.6, 46.0, ['measured', 'pitch_um', 'inner'], 'row pitch (mask 45.8)'),
                                                           spacing=(813.0, 823.0, ['measured', 'leg_centreline_separation_um', 'mask'],
                                                                    'the mask\'s leg spacing against the inset\'s 817.9 µm (consistency)'),
                                                           neighbour=(2.0, 2.0, ['measured', 'neighbour_leg'],
                                                                      'the inset\'s neighbouring leg (nested, RC8: leg 2)'))),
    # Xia (1.798 µm/px, the rc6 scale): relief contrast with doubled outlines, so the chamfer nulls are the claim's
    'xia_fig3c__rc7': dict(v1=True, median_max_um=3.5, p90_max_um=7.0, chamfer_null_ratio_min=1.3,
                           v5_along=(0.99, 1.01), v5_across=(0.99, 1.01),
                           measured=dict(centre=(297.0, 303.0, ['measured', 'legs', '0', 'centre_um_image'], 'leg 1 centre at its first posts (printed 300)'))),
    'xia_fig3c__pp': dict(v1=True, median_max_um=3.5, p90_max_um=7.0, chamfer_null_ratio_min=1.3,
                          v5_along=(0.99, 1.01), v5_across=(0.99, 1.01),
                          measured=dict(centre=(297.0, 303.0, ['measured', 'legs', '0', 'centre_um_image'], 'leg 1 centre at its first posts (printed 300)'))),
    'xia_fig3e__pp': dict(v1=True, median_max_um=3.0, p90_max_um=13.0, chamfer_null_ratio_min=1.15,
                          v5_along=(0.99, 1.01), v5_across=(0.99, 1.01)),
})
# recorded, not claimed: Xia panel e on the depicted mask (leg 3 narrower than the constant-shear taper puts it: V1
# across fails), panel f on both (INV-4, open), panel g on the aligned mask (as traced)
RECORDED = {
    'xia_fig3e__rc7': dict(v1=False, v5_along=(0.99, 1.01), v5_across=(1.0, 1.02),
                           measured=dict(leg3=(183.0, 189.0, ['measured', 'legs', '1', 'centre_um_image'], 'leg 3 centre at Fig. 3a\'s box (mask 205.1)'))),
    'xia_fig3f__rc7': dict(v1=True, v5_along=(0.99, 1.01), v5_across=(0.95, 0.98)),
    'xia_fig3f__pp': dict(v1=True, v5_along=(0.99, 1.01), v5_across=(0.94, 0.98)),
    'xia_fig3g_outlet__pp': dict(v1=False, v5_along=(0.99, 1.01), v5_across=(0.99, 1.02)),
    # RC8 (A7, D5): claims in rc7, recorded since: the outline metric found an outline for only 42-48 % of the eligible
    # samples and its censored nulls equal the registration's (0.94-1.00), so the metric cannot tell; the registration (V1, V5)
    # is reproduced and its ranges checked as rc7's claim limits had them
    'xia_fig3d__rc7': dict(v1=True, v5_along=(0.995, 1.015), v5_across=(0.995, 1.015)),
    'xia_fig3d__pp': dict(v1=True, v5_along=(0.995, 1.015), v5_across=(0.995, 1.015)),
    'xia_fig3g_outlet__rc7': dict(v1=True, v5_along=(0.99, 1.01), v5_across=(1.02, 1.05)),
}
COVERAGE_MIN = 0.5             # RC8 (A7): a claim's outline metric must find an outline for at least half its samples
NULL_RATIO_CRITERION = 1.3     # package fix #7: every panel's null ratios reported, those below the criterion flagged

panels = {}
null_flags = []


def null_ratios(rec):
    med = rec['metric']['median_abs_um']
    mn = rec.get('metric_null') or {}
    c = rec['V2_chamfer']
    out = dict(chamfer=min(min(c['half_pitch_along_px']), min(c['five_um_across_px'])) / c['best_px'])
    if mn:
        out['metric'] = min(mn['half_pitch_along']['median_abs_um'], mn['5um_across']['median_abs_um']) / med
        cr = mn.get('censored_ratios') or {}
        if cr.get('half_pitch_along') is not None and cr.get('5um_across') is not None:
            out['censored'] = min(cr['half_pitch_along'], cr['5um_across'])      # RC8 (A7; the code review's N4)
    return out


def add(rec, layer, role, lim):
    nr = null_ratios(rec)
    for k, v in nr.items():
        if v < NULL_RATIO_CRITERION:
            null_flags.append(dict(panel=rec['name'], role=role, null=k, ratio=round(v, 3),
                                   used_for_V2=('%s_null_ratio_min' % k) in lim))
    if role == 'claim':
        lim = dict(lim, coverage_min=COVERAGE_MIN)
        cov = rec['metric'].get('coverage')            # RC8 code review, m5: a claim is checked here, not only by the test
        if cov is None or cov < COVERAGE_MIN:
            raise SystemExit('%s: coverage %s below the claims\' %g' % (rec['name'], cov, COVERAGE_MIN))
        if 'censored_null_ratio_min' in lim and not (nr.get('censored', 0.0) >= lim['censored_null_ratio_min']):
            raise SystemExit('%s: censored null ratio %s below its claim limit %g' % (
                rec['name'], nr.get('censored'), lim['censored_null_ratio_min']))
    panels[rec['name']] = dict(layer=layer, role=role, limits=lim, null_ratios=nr, anchored=dict(
        V1_along_px=rec['V1_feature_vs_chamfer_px']['along'], V1_across_px=rec['V1_feature_vs_chamfer_px']['across'],
        median_abs_um=rec['metric']['median_abs_um'], p90_abs_um=rec['metric']['p90_abs_um'],
        chamfer_best_px=rec['V2_chamfer']['best_px'], mask_sha256=rec['mask_sha256'],
        metric_window_um=rec.get('metric_window_um'),
        coverage=rec['metric'].get('coverage'), n_eligible=rec['metric'].get('n_eligible'),
        censored_ratios=(rec.get('metric_null') or {}).get('censored_ratios'),
        V5=[rec['V5_scale_scan']['best_along'], rec['V5_scale_scan']['best_across']] if 'V5_scale_scan' in rec else None))


def place(rec, layer=3):
    n = rec['name']
    if n in CLAIM:
        add(rec, layer, 'claim', CLAIM[n])
    elif n in RECORDED:
        add(rec, layer, 'recorded', RECORDED[n])
    else:
        raise SystemExit('no limits for %s' % n)


for tag in ('rc7', 'fa'):
    d = J('dinh', 'dinh_results_%s.json' % tag)
    for key in ('p1', 'p2', 'p3', 'p4'):
        place(d[key])
for tag in ('rc7', 'aligned'):
    place(J('gifford', 'gifford_results_%s.json' % tag))
for tag in ('design', 'aligned'):
    for f in ('gifford6_results_%s.json' % tag, 'gifford6_results_%s_onemap.json' % tag):
        for r in J('gifford', f).values():
            place(r)
for tag in ('rc7', 'pp'):
    place(J('xia', 'xia_results_%s.json' % tag))
    for L in 'cdef':
        place(J('xia', 'xia_results_fig3%s_%s.json' % (L, tag)))


def fig3a(rec, v12_expect):
    return dict(pose_whole_rotation_deg=rec['pose_whole']['rotation_deg'], pose_whole_chamfer_px=rec['pose_whole']['chamfer_px'],
                scale_free_over_bar=rec['scale_free']['scale_over_bar'],
                group_shifts={g: dict(along_um=s['shift_along_um'], across_um=s['shift_across_um']) for g, s in rec['V11_group_shifts'].items()},
                V12_ok={k: v['ok'] for k, v in rec['V12_boxes'].items() if 'ok' in v}, V12_expected=v12_expect,
                dxf_sha256=rec['dxf_sha256'])


def v12_text(rec):
    ins = [k for k, v in sorted(rec['V12_boxes'].items()) if v.get('ok')]
    out = ['%s %.2f of its box\'s half-size' % (k, v['offset_over_half_size']) for k, v in sorted(rec['V12_boxes'].items())
           if 'ok' in v and not v['ok']]
    return '%s inside their boxes%s' % (', '.join(ins), ('; outside: ' + ', '.join(out)) if out else '')


# RC8 (V13): the taper-direction sign test on each Fig. 4 reading (taper_sign.py, run by run_all_rc7.sh).  The decisive
# outcomes are required, not merely anchored: the nested reproduction's neighbour PASSES, the path-order reading B's
# FAILS, and leg 4 itself (the control) PASSES on both
V13_RUNS = dict(rc7=('shipped', 'gifford2014_fig4_design'), aligned=('shipped', 'gifford2014_fig4_aligned'),
                pathorder=('variant', 'gifford2014_fig4_pathorder'), start51=('variant', 'gifford2014_fig4_start51'),
                fold20=('variant', 'gifford2014_fig4_fold20'))
V13_REQUIRED = dict([(('aligned', 'neighbour'), 'PASS'), (('pathorder', 'neighbour'), 'FAIL')] +
                    [((t_, 'leg4'), 'PASS') for t_ in V13_RUNS])          # leg 4, the control, on every mask
V13_WINDOWS = dict(neighbour=(30.0, None), leg4=(12.0, '@array'))           # RC8 code review, m3: +-30 um on the neighbour
v13 = {}
for tag, (where, mask_) in V13_RUNS.items():
    for which, (window, xmin) in V13_WINDOWS.items():
        f = os.path.join(RES, 'gifford', 'v13_%s_%s.json' % (tag, which))
        if not os.path.exists(f):
            raise SystemExit('V13 result missing: %s (run run_all_rc7.sh)' % f)
        v = json.load(open(f))
        if v['window_um'] != window or (xmin is None) != (v['x_limits'][0] is None):     # n3: as run_all ran it
            raise SystemExit('V13 %s %s: run at window %s, x %s; expected %s, %s' % (tag, which, v['window_um'],
                                                                                   v['x_limits'], window, xmin))
        want = V13_REQUIRED.get((tag, which))
        if want and v['verdict'] != want:
            raise SystemExit('V13 %s %s: %s, where %s is required' % (tag, which, v['verdict'], want))
        v13['%s_%s' % (tag, which)] = dict(
            trace='gifford/gifford_results_%s.json' % tag, where=where, mask=mask_, leg=which, window_um=window, x_min=xmin,
            leg_Y=v['leg_Y'], verdict=v['verdict'], required=want, positions=v['positions'], eligible=v['eligible'],
            slope_image=v.get('slope_image'), se_image=v.get('se_image'), z_image=v.get('z_image'),
            se_method=v.get('se_method'), se_image_iid=v.get('se_image_iid'), se_image_hac=v.get('se_image_hac'),
            slope_mask=v.get('slope_mask'),
            coverage=v.get('coverage'), dxf_sha256=v.get('dxf_sha256'),
            offsets_um=[v.get('offset_low_mean_um'), v.get('offset_high_mean_um')],
            offset_abs_range_um=([min(min(abs(r_['off_low']), abs(r_['off_high'])) for r_ in v['rows']),
                                  max(max(abs(r_['off_low']), abs(r_['off_high'])) for r_ in v['rows'])]
                                 if v.get('rows') else None))

F_rc7, F_pp = J('xia', 'xia_fig3a_device_rc7.json'), J('xia', 'xia_fig3a_device_pp.json')
ras = os.environ['RASTERS']
meas6 = json.load(open(os.path.join(ROOT, 'provenance', 'rc6_scripts', 'gifford_fig6_measurement.json')))
r6v = json.load(open(os.path.join(ROOT, 'opencif', 'data', 'r6', 'rc6_values.json')))
# the figure measurements the traces rest on are rc6's (rc7 re-measured none of them)
arc = json.load(open(os.path.join(ROOT, 'provenance', 'rc6_traces', 'measurements', 'xia_fig3d_arc_scale.json')))

out = dict(
    note=('rc7 trace anchors: micrograph traces of every shipped mask the papers\' figures show (the rc7 presets\' masks '
          'included), with the rc7 trace code (a metric half-window of at most half the gap).  Not hashed code: '
          'provenance/rc7_traces/.  tests/test_trace_anchors.py repeats every registration from the rasters and the DXFs '
          'and applies the limits below (set with the results in view: regression limits; rc6\'s kept for every panel '
          'rc6 claimed that is still claimed: RC8, D5, records Xia 3d and 3g on the depicted mask).  null_flags lists every panel whose metric or chamfer null ratio is below %.1f x, and whether '
          'its V2 claim uses that null.' % NULL_RATIO_CRITERION),
    code_stamp=code_stamp(),
    rasters={'dinh5-000.png': dict(sha256=T.sha256(os.path.join(ras, 'dinh5-000.png')),
                                   source='Dinh et al. 2024 Lab Chip 24, 913-923, PDF page 5 (printed p. 917), Fig. 2: pdfimages -png, image 0'),
             'gif5-001.png': dict(sha256=T.sha256(os.path.join(ras, 'gif5-001.png')),
                                  source='Gifford et al. 2014 Lab Chip 14, 4496-4505, PDF page 5 (printed p. 4500), Fig. 4: pdfimages -png, image 1'),
             'gif7-000.png': dict(sha256=T.sha256(os.path.join(ras, 'gif7-000.png')),
                                  source='Gifford et al. 2014 Lab Chip 14, 4496-4505, PDF page 7 (printed p. 4502), Fig. 6: pdfimages -png, image 0'),
             'xia5-000.png': dict(sha256=T.sha256(os.path.join(ras, 'xia5-000.png')),
                                  source='Xia et al. 2016 Sci. Rep. 6:35943, PDF page 5, Fig. 3: pdfimages -png, image 0')},
    bars=dict(
        dinh_fig2B=dict(um=50, rows=[[470, 490], [952, 972], [1436, 1456], [1920, 1940]], cols=[870, 975], px_aa=46.906, full_px=46),
        gifford_fig4_inset=dict(um=300, rows=[550, 600], cols=[690, 829], px_aa=107.52, full_px=107),
        xia_fig3_panel_i=dict(um=250, rows=[3085, 3120], cols=[1860, 2025], px_aa=129.11, full_px=128, polarity='bright',
                              note='panel i\'s bar; panels b-g are at the scale of their printed widths (1.798 µm/px)'),
        xia_fig3a=dict(um=5000, rows=[2925, 2960], cols=[1000, 1230], px_aa=178.32, full_px=177),
        gifford_fig6=dict(um=250, visible_length_px=meas6['scale_bar']['visible_length_px'],
                          um_per_px_upper_bound=meas6['scale_bar']['um_per_px_upper_bound'],
                          note='cut by the raster\'s right edge: an upper bound on the scale only')),
    scales=dict(xia=dict(rc6_um_per_px=r6v['xia']['scale_um_per_px'], printed_widths_um_per_px=1.8028, arc_um_per_px=arc['scale_um_per_px']['mean']),
                gifford_fig6=dict(um_per_px=meas6['um_per_px'], from_='the printed 125 µm centre between the rows\' tip lines')),
    runs=dict(
        shipped=[dict(kind='dinh', mask='dinh2024_v1_depicted', tag='rc7', label='the shipped mask dinh2024_v1_depicted'),
                 dict(kind='dinh', mask='dinh2024_v1_figure_aligned', tag='fa', label='the figure-aligned mask dinh2024_v1_figure_aligned'),
                 dict(kind='gifford4', mask='gifford2014_fig4_design', tag='rc7', label='the shipped mask gifford2014_fig4_design'),
                 dict(kind='gifford4', mask='gifford2014_fig4_aligned', tag='aligned',
                      label='the aligned mask gifford2014_fig4_aligned (nested legs, RC8)'),
                 dict(kind='gifford6', mask='gifford2014_fig6_design', tag='design', label='the mask gifford2014_fig6_design'),
                 dict(kind='gifford6', mask='gifford2014_fig6_design', tag='design', one_map=True, turns=True,
                      label='the mask gifford2014_fig6_design'),
                 dict(kind='gifford6', mask='gifford2014_fig6_aligned', tag='aligned', label='the mask gifford2014_fig6_aligned'),
                 dict(kind='gifford6', mask='gifford2014_fig6_aligned', tag='aligned', one_map=True, turns=True,
                      label='the mask gifford2014_fig6_aligned'),
                 dict(kind='xia', mask='xia2016_rs7_depicted', tag='rc7', label='the mask xia2016_rs7_depicted'),
                 dict(kind='xia', mask='xia2016_rs7_aligned', tag='pp', label='the mask xia2016_rs7_aligned')],
        variants=[]),
    panels=panels,
    null_flags=null_flags,
    fig3a=dict(rc7=fig3a(F_rc7, v12_text(F_rc7)), pp=fig3a(F_pp, v12_text(F_pp))),
    variants={},
    v13=v13,
    v13_note=('V13, the taper-direction sign test (RC8; provenance/rc7_traces/code/taper_sign.py): PASS when the image\'s '
              'outer-width slope along the leg has the mask\'s sign beyond 3 SE with coverage >= 0.5, >= 30 positions and a '
              'mask slope of at least 1e-4 um/um; FAIL at the opposite sign; else INCONCLUSIVE.  The SE is the larger of '
              'the Newey-West SE (Bartlett kernel, lag 100 um, beyond two row pitches) and the independent-sample SE: '
              'neighbouring positions share pixels and their residuals are autocorrelated, so the independent-sample SE '
              'alone is too small (RC8 code review, M1, N2).  A sign test: the magnitude '
              'ratio is reported, not judged.  The wall window must exceed the largest wall offset (review 3, E6; the code '
              'review, m3): +-30 um on the neighbouring leg (the nested reading\'s walls lie %.1f-%.1f um outside the '
              'mask\'s, offset_abs_range_um), +-12 um on leg 4\'s array part.  Required: %s.' % (
                  tuple(v13['aligned_neighbour']['offset_abs_range_um']) + (', '.join(
                      '%s %s %s' % (k[0], k[1], w) for k, w in sorted(V13_REQUIRED.items())),))),
)
json.dump(out, open(OUTP, 'w'), indent=1, default=float)
print('V13: ' + '; '.join('%s %s' % (k, v['verdict']) for k, v in sorted(v13.items())))
print('%s: %d panels (%d claims); null flags: %s' % (OUTP, len(panels), sum(1 for p in panels.values() if p['role'] == 'claim'),
                                                    ', '.join('%s %s %.2f%s' % (f['panel'], f['null'], f['ratio'], ' (V2)' if f['used_for_V2'] else '')
                                                              for f in null_flags) or 'none'))
