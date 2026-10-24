"""Collect the rc6 trace records into provenance/rc6_traces/trace_anchors.json: the inputs (raster hashes, the bars read
with their anti-aliased ends, the scales), the runs tests/test_trace_anchors.py repeats, each panel's anchored
registration and metric, the acceptance limits of the claims and the ranges of the recorded panels, Xia's Fig. 3a
whole-device registrations and the variant presets' declared fields.

The limits were set with the results in view: they are regression limits for rc6, not a blind test.  What is blind is
each panel's scale (a bar, a printed width or the measurement), its rotation (the image) and the agreement of two
independent registrations (V1).

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


# ---- the claims' limits and the recorded panels' ranges (V5 ranges only where the scan ran)
# Dinh (1.066 µm/px, the 50 µm bars): every panel a claim on the metric's nulls
CLAIM = {
    'dinh_fig2B_p1_head__rc6': dict(v1=True, median_max_um=1.6, p90_max_um=3.2, metric_null_ratio_min=1.3,
                                    v5_along=(0.97, 1.0), v5_across=(0.975, 1.0)),
    'dinh_fig2B_p2_transition__rc6': dict(v1=True, median_max_um=1.0, p90_max_um=3.0, metric_null_ratio_min=2.0,
                                          v5_along=(0.98, 1.0), v5_across=(0.995, 1.015)),
    'dinh_fig2B_p3_linear__rc6': dict(v1=True, median_max_um=0.6, p90_max_um=1.2, metric_null_ratio_min=2.0,
                                      v5_along=(0.99, 1.01), v5_across=(0.99, 1.01),
                                      measured=dict(centre=(68.0, 69.5, ['measured', 'image_widths_um', 'centre'], 'centre width (mask 68.1 at its row)'))),
    'dinh_fig2B_p4_outlet__rc6': dict(v1=True, median_max_um=0.7, p90_max_um=2.0, metric_null_ratio_min=2.0,
                                      v5_along=(0.99, 1.01), v5_across=(0.99, 1.01),
                                      measured=dict(centre=(51.5, 53.0, ['measured', 'image_widths_um', 'centre'], 'centre width (printed 52)'),
                                                    side=(43.5, 44.5, ['measured', 'image_widths_um', 'side_L'], 'left side width (printed 44)'),
                                                    pitch=(44.9, 45.3, ['measured', 'pitch_um'], 'row pitch (mask 45.1)'))),
    # Gifford's Fig. 4 inset (2.79 µm/px, the 300 µm bar): relief contrast, outlines read long along the flow
    'gifford_fig4_inset_leg4__rc6': dict(v1=True, median_max_um=1.5, p90_max_um=20.0, metric_null_ratio_min=2.5,
                                         v5_along=(0.99, 1.01), v5_across=(0.99, 1.01),
                                         measured=dict(pitch=(45.6, 46.0, ['measured', 'pitch_um', 'inner'], 'row pitch (mask 45.8)'),
                                                       spacing=(813.0, 823.0, ['measured', 'leg_centreline_separation_um', 'mask'],
                                                                'the mask\'s leg spacing: the inset\'s 817.9 µm (rc6 preset; consistency)'))),
    # Gifford's Fig. 6 (4.52 µm/px from the printed centre; each leg on its own map): the side channels read 5-9 µm wider
    # than drawn (reported), so the metric's 5 µm null is close to the median; the chamfer nulls are the claim's
    # (V1 and the chamfer's across nulls; along the flow the pairing rests on END_CHECK, which the test checks separately)
    **{'gifford_fig6_leg%d__rc6' % j: dict(v1=True, median_max_um=2.0, p90_max_um=10.0, chamfer_null_ratio_min=1.15,
                                           v5_along=(0.995, 1.005), v5_across=(1.0, 1.05),
                                           v5_note='along: consistency, the pitch was measured on this image; across: the side channels read wider')
       for j in (7, 5, 3, 1)},
    # the preset takes the micrograph's leg axes (rc6, after the 2026-10-17 review): the four legs and the turns of
    # legs 5, 3, 1 on ONE map
    'gifford_fig6_four_legs__rc6': dict(v1=True, median_max_um=2.0, p90_max_um=8.0, chamfer_null_ratio_min=1.15,
                                                       v5_along=(0.995, 1.005), v5_across=(0.995, 1.005),
                                                       v5_note='consistency: the pitch and the leg axes were read on this image',
                                                       measured=dict(
        leg7=(-1.0, 1.0, ['measured', 'leg7_own_map_minus_common_um', 'across'],
              'leg 7 (not in the fit) across, from its own end\'s map, µm (consistency: the axes were read on this image)'),
        # the turns: the leg axes were read on these turn walls, so distance, radii and the centre across the rows are
        # consistency (they check that the generator draws the axes it is given); the centre along the flow is not set by
        # anything read on the turns (the generator places the turns from the legs' ends): the one independent turn result
        turn_med=(0.0, 2.0, ['measured', 'turns', 'nearest_mask_edge_um', 'median'],
                  'the turns: median distance from the image\'s turn walls to the mask\'s nearest edge, µm (consistency; the first, uniform form: 23.6)'),
        turn_p90=(0.0, 5.0, ['measured', 'turns', 'nearest_mask_edge_um', 'p90'],
                  'the turns: 90th percentile of that distance, µm (consistency; the first form: 59.7)'),
        turn_r=(-2.0, 2.0, ['measured', 'turns', 'radii_about_image_centre', 'mean_um'],
                'the turns: mean radius mask - image about the image\'s turn centre, µm (consistency)'),
        turn_c=(-20.0, 20.0, ['measured', 'turns', 'centre_offset_mask_minus_image_um', 'v'],
                'the turns: the mask\'s own turn centre across from the image\'s, µm (consistency; the first form: -1454)'),
        turn_u=(-10.0, 10.0, ['measured', 'turns', 'centre_offset_mask_minus_image_um', 'u'],
                'the turns: the mask\'s own turn centre along the flow from the image\'s, µm (+ upstream: -5.7 is 5.7 µm '
                'downstream; the image centre\'s SD 0.5 µm): the one independent turn result'))),
    # Xia (1.798 µm/px, the rc6 scale): relief contrast with doubled outlines, so the chamfer nulls are the claim's
    'xia_fig3c__rc6': dict(v1=True, median_max_um=3.5, p90_max_um=7.0, chamfer_null_ratio_min=1.3,
                           v5_along=(0.99, 1.01), v5_across=(0.99, 1.01),
                           measured=dict(centre=(297.0, 303.0, ['measured', 'legs', '0', 'centre_um_image'], 'leg 1 centre at its first posts (printed 300)'))),
    'xia_fig3d__rc6': dict(v1=True, median_max_um=3.0, p90_max_um=10.0, chamfer_null_ratio_min=1.3,
                           v5_along=(0.995, 1.015), v5_across=(0.995, 1.015)),
    'xia_fig3g_outlet__rc6': dict(v1=True, median_max_um=3.0, p90_max_um=13.0, chamfer_null_ratio_min=1.04,
                                  v5_along=(0.99, 1.01), v5_across=(1.02, 1.05)),
}
# recorded, not claimed: Xia panel e (leg 3 narrower than the constant-shear taper puts it: V1 across fails) and panel f
# (at its own pitch scale; V5 across 0.96)
RECORDED = {
    'xia_fig3e__rc6': dict(v1=False, v5_along=(0.99, 1.01), v5_across=(1.0, 1.02),
                           measured=dict(leg3=(183.0, 189.0, ['measured', 'legs', '1', 'centre_um_image'], 'leg 3 centre at Fig. 3a\'s box (mask 205.1)'))),
    'xia_fig3f__rc6': dict(v1=True, v5_along=(0.99, 1.01), v5_across=(0.95, 0.98)),
}
VARIANT_LIMITS = {
}

panels = {}


def add(rec, layer, role, lim):
    panels[rec['name']] = dict(layer=layer, role=role, limits=lim, anchored=dict(
        V1_along_px=rec['V1_feature_vs_chamfer_px']['along'], V1_across_px=rec['V1_feature_vs_chamfer_px']['across'],
        median_abs_um=rec['metric']['median_abs_um'], p90_abs_um=rec['metric']['p90_abs_um'],
        chamfer_best_px=rec['V2_chamfer']['best_px'], mask_sha256=rec['mask_sha256'],
        V5=[rec['V5_scale_scan']['best_along'], rec['V5_scale_scan']['best_across']] if 'V5_scale_scan' in rec else None))


def place(rec, layer, table, other):
    n = rec['name']
    if n in table:
        add(rec, layer, 'claim', table[n])
    elif n in other:
        add(rec, layer, 'recorded', other[n])
    else:
        raise SystemExit('no limits for %s' % n)


for key in ('p1', 'p2', 'p3', 'p4'):
    place(J('dinh', 'dinh_results_rc6.json')[key], 3, CLAIM, RECORDED)
place(J('gifford', 'gifford_results_rc6.json'), 3, CLAIM, RECORDED)
for f in ('gifford6_results_rc6.json', 'gifford6_results_rc6_onemap.json'):
    for r in J('gifford', f).values():
        place(r, 3, CLAIM, RECORDED)
place(J('xia', 'xia_results_rc6.json'), 3, CLAIM, RECORDED)
for L in 'cdef':
    place(J('xia', 'xia_results_fig3%s_rc6.json' % L), 3, CLAIM, RECORDED)
# the Xia variants: recorded (their claim is the whole-device registration on Fig. 3a)
for f in (['xia_results_dfill.json'] + ['xia_results_fig3%s_dfill.json' % L for L in 'cdef'] +
          ['xia_results_fig3.json'] + ['xia_results_fig3%s_fig3.json' % L for L in 'cdef']):
    r = J('xia', f)
    v5 = r.get('V5_scale_scan') or {}
    add(r, 5, 'recorded', dict(v1=bool(r['V1_feature_vs_chamfer_px']['agree'])))


def fig3a(rec, v12_expect):
    return dict(pose_whole_rotation_deg=rec['pose_whole']['rotation_deg'], pose_whole_chamfer_px=rec['pose_whole']['chamfer_px'],
                scale_free_over_bar=rec['scale_free']['scale_over_bar'],
                group_shifts={g: dict(along_um=s['shift_along_um'], across_um=s['shift_across_um']) for g, s in rec['V11_group_shifts'].items()},
                V12_ok={k: v['ok'] for k, v in rec['V12_boxes'].items() if 'ok' in v}, V12_expected=v12_expect,
                dxf_sha256=rec['dxf_sha256'])


F_rc6, F_df, F_f3 = J('xia', 'xia_fig3a_device_rc6.json'), J('xia', 'xia_fig3a_device_dfill.json'), J('xia', 'xia_fig3a_device_fig3.json')
FIG3_V12 = 'as anchored (V12_ok); see the rc6_traces README'
ras = os.environ['RASTERS']
meas6 = json.load(open(os.path.join(ROOT, 'provenance', 'rc6_scripts', 'gifford_fig6_measurement.json')))
r6v = json.load(open(os.path.join(ROOT, 'opencif', 'data', 'r6', 'rc6_values.json')))
arc = json.load(open(os.path.join(os.path.dirname(RES.rstrip('/')), 'measurements', 'xia_fig3d_arc_scale.json')))


def rec_rows(p):
    d = json.load(open(p))['design']
    return d['N'], d['NM']


VAR = os.environ.get('OPENCIF_TRACE_VARIANTS', os.path.join(ROOT, 'provenance', 'rc6_traces', 'variants', 'runs'))
_fit = J('xia', 'xia_fig3_fit.json')['result']
FIG3_DECLARED = {'targets.cutoff_um': _fit['D_um'], 'targets.filtrate_fraction': _fit['F']}
variants = {}
for name, base, declared in (('xia2016_rs7_depicted_dfill', 'xia2016_rs7_depicted', {'targets.cutoff_um': 6.531498914336713}),
                             ('xia2016_rs7_depicted_fig3', 'xia2016_rs7_depicted', FIG3_DECLARED)):
    N, NM = rec_rows(os.path.join(VAR, name, name + '.json'))
    variants[name] = dict(base=base, declared=declared, N=N, NM=NM)

out = dict(
    note=('rc6 trace anchors (2026-10-17): micrograph traces of the shipped rc6 masks and two variants on the papers\' figures.  '
          'Not hashed code: provenance/rc6_traces/ (README.md there).  tests/test_trace_anchors.py repeats every registration '
          'from the rasters and the DXFs and applies the limits below (set with the results in view: regression limits).'),
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
        shipped=[dict(kind='dinh', mask='dinh2024_v1_depicted', tag='rc6', label='the shipped mask dinh2024_v1_depicted'),
                 dict(kind='gifford4', mask='gifford2014_fig4_design', tag='rc6', label='the shipped mask gifford2014_fig4_design'),
                 dict(kind='gifford6', mask='gifford2014_fig6_design', tag='rc6', label='the shipped mask gifford2014_fig6_design'),
                 dict(kind='gifford6', mask='gifford2014_fig6_design', tag='rc6', one_map=True, turns=True,
                      label='the shipped mask gifford2014_fig6_design (the micrograph\'s leg axes)'),
                 dict(kind='xia', mask='xia2016_rs7_depicted', tag='rc6', label='the shipped mask xia2016_rs7_depicted')],
        variants=[dict(kind='xia', mask='xia2016_rs7_depicted_dfill', tag='dfill',
                       label='the variant mask xia2016_rs7_depicted_dfill (D 6.5315 um)'),
                  dict(kind='xia', mask='xia2016_rs7_depicted_fig3', tag='fig3',
                       label='the variant mask xia2016_rs7_depicted_fig3 (F 0.94072, D 6.5705 um, read off Figs. 3g and 3a)')]),
    panels=panels,
    fig3a=dict(rc6=fig3a(F_rc6, 'c, d, e, f inside their boxes; g outside (the drawing\'s outlet lies 5.4 mm beyond the mask\'s)'),
               dfill=fig3a(F_df, 'c, d, e, f inside their boxes; g 1.08 of its box\'s half-size (0.46 mm of its 0.65 mm offset is '
                                 'panel g\'s own offset from its box, which the shipped mask shows)'),
               fig3=fig3a(F_f3, FIG3_V12)),
    variants=variants,
)
json.dump(out, open(OUTP, 'w'), indent=1, default=float)
print('%s: %d panels (%d claims)' % (OUTP, len(panels), sum(1 for p in panels.values() if p['role'] == 'claim')))
