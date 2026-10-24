"""The results table of provenance/rc6_traces/README.md, from the trace records.
usage: python3 summary_table.py RESULTS_DIR OUT.md"""
import json
import os
import sys

RES, OUT = sys.argv[1:3]


def J(*p):
    return json.load(open(os.path.join(RES, *p)))


rows = []


def row(label, mask, r, role):
    v1 = r['V1_feature_vs_chamfer_px']
    m, n, c = r['metric'], r['metric_null'], r['V2_chamfer']
    v5 = r.get('V5_scale_scan')
    rows.append('| %s | %s | %.4f | %.2f / %.2f | %.2f / %.2f | %.2f / %.2f | %.2f / %.2f / %.2f | %s | %s |' % (
        label, mask, 1.0 / r['px_per_um'], v1['along'], v1['across'], m['median_abs_um'], m['p90_abs_um'],
        n['half_pitch_along']['median_abs_um'], n['5um_across']['median_abs_um'],
        c['best_px'], min(c['half_pitch_along_px']), min(c['five_um_across_px']),
        (('%.3f / %.3f' % (v5['best_along'], v5['best_across'])) + (' (grid edge)' if v5.get('on_grid_edge') else '')) if v5 else '—',
        role))


d = J('dinh', 'dinh_results_rc6.json')
for k, lab in (('p1', 'Dinh 2B p1 head'), ('p2', 'Dinh 2B p2 transition'), ('p3', 'Dinh 2B p3 linear'), ('p4', 'Dinh 2B p4 outlet')):
    row(lab, 'shipped', d[k], 'claim')
row('Gifford Fig. 4 inset, leg 4', 'shipped (design)', J('gifford', 'gifford_results_rc6.json'), 'claim')
g6 = J('gifford', 'gifford6_results_rc6.json')
for j in (7, 5, 3, 1):
    row('Gifford Fig. 6, leg %d' % j, 'shipped (design)', g6['leg%d' % j], 'claim')
row('Gifford Fig. 6, legs 7-5-3-1 and their turns on one map', 'shipped (design)',
    J('gifford', 'gifford6_results_rc6_onemap.json')['four_legs'], 'claim')
for f, lab, role in (('xia_results_fig3c_rc6.json', 'Xia 3c leg 1, first posts', 'claim'),
                     ('xia_results_fig3d_rc6.json', 'Xia 3d legs 5, 3 below the loops', 'claim'),
                     ('xia_results_fig3e_rc6.json', 'Xia 3e legs 5, 3 half way', 'recorded'),
                     ('xia_results_fig3f_rc6.json', 'Xia 3f leg 7 (own scale)', 'recorded'),
                     ('xia_results_rc6.json', 'Xia 3g outlet', 'claim')):
    row(lab, 'shipped (depicted)', J('xia', f), role)
for tag, mlab in (('dfill', 'variant dfill'), ('fig3', 'variant fig3 (F, D off Fig. 3)')):
    for f, lab in (('xia_results_fig3c_%s.json', 'Xia 3c'), ('xia_results_fig3d_%s.json', 'Xia 3d'),
                   ('xia_results_fig3e_%s.json', 'Xia 3e'), ('xia_results_fig3f_%s.json', 'Xia 3f'), ('xia_results_%s.json', 'Xia 3g')):
        p = os.path.join(RES, 'xia', f % tag)
        if os.path.exists(p):
            row(lab, mlab, json.load(open(p)), 'recorded')
head = ('| panel | mask | µm/px | V1 along / across px | median / p90 µm | metric null ½ pitch / 5 µm | chamfer best / ½ pitch / '
        '5 µm px | V5 along / across | role |\n|---|---|---|---|---|---|---|---|---|')
F = []
for tag in ('rc6', 'dfill', 'fig3'):
    p = os.path.join(RES, 'xia', 'xia_fig3a_device_%s.json' % tag)
    if os.path.exists(p):
        f = json.load(open(p))
        sh = f['V11_group_shifts']
        F.append('| Xia Fig. 3a, %s | %.3f° / %.3f px | %+.0f / %+.0f µm | %+.0f / %+.0f µm | %s |' % (
            {'rc6': 'shipped', 'dfill': 'variant dfill', 'fig3': 'variant fig3'}[tag], f['pose_whole']['rotation_deg'], f['pose_whole']['chamfer_px'],
            sh['loops at the bottom']['shift_along_um'], sh['loops at the bottom']['shift_across_um'],
            sh['outlet ring and ports']['shift_along_um'], sh['outlet ring and ports']['shift_across_um'],
            ', '.join('%s %s' % (k, 'ok' if v['ok'] else 'OUT (%.1f)' % v['offset_over_half_size'])
                      for k, v in f['V12_boxes'].items() if 'ok' in v)))
fh = ('| whole device | pose (whole device) rotation / chamfer | far loops along / across | outlet ring along / across | V12 panels in their boxes |\n'
      '|---|---|---|---|---|')
# Gifford Fig. 6, the turns (panels_gifford6.turn_compare): the first variant (uniform 588 µm legs) and the current one
T_ = []
for f, lab in (('gifford6_results_fig6legs588.json', 'uniform legs 588 µm apart (first form, variants/history588)'),
               ('gifford6_results_rc6_onemap.json', 'the preset, with the measured leg axes (shipped)')):
    p = os.path.join(RES, 'gifford', f)
    if os.path.exists(p):
        t = json.load(open(p))['four_legs']['measured'].get('turns')
        if t:
            co = t['centre_offset_mask_minus_image_um']
            conc = (co['u'] ** 2 + co['v'] ** 2) ** 0.5 <= 50.0        # the mask's turns concentric with the image's
            # the flow runs along -u (panels_gifford6): -u is the distance downstream
            T_.append('| %s | %d / %d | %.1f / %.1f / %.0f µm | %+.1f / %+.0f µm | %s |' % (
                lab, t['image']['walls'], t['mask']['walls'], t['nearest_mask_edge_um']['median'], t['nearest_mask_edge_um']['p90'],
                t['nearest_mask_edge_um']['max'], -co['u'], co['v'],
                ('%+.1f ± %.1f µm' % (t['radii_about_image_centre']['mean_um'], t['radii_about_image_centre']['sd_um'])) if conc
                else '— (not concentric with the image\'s turns)'))
th = ('| Gifford Fig. 6 turns, mask | walls image / mask | image wall to nearest mask edge: median / p90 / max | turn centre mask - '
      'image, downstream / across | radii mask - image (about the image\'s centre) |\n|---|---|---|---|---|')
tn = ('\nThe shipped preset\'s leg axes are read on these turns, so for it the nearest-edge distances, the centre across the rows and '
      'the radii agree by construction (consistency); the centre\'s downstream offset is the one independent turn result.  '
      '"(grid edge)": the V5 scan\'s best scale lies on the edge of its grid (a bound, not a minimum).\n')
open(OUT, 'w').write(head + '\n' + '\n'.join(rows) + '\n\n' + fh + '\n' + '\n'.join(F) + '\n' +
                     (('\n' + th + '\n' + '\n'.join(T_) + '\n') if T_ else '') + tn)
print(open(OUT).read())
