"""EVIDENCE.md: for every paper figure the rc6 traces cover, where its evidence is -- the shipped preset, the zoomed trace
(each panel at its own scale), the whole-device figure (every polygon of the DXF with the micrographs placed by their
stored maps), the registration's numbers -- and what that evidence tests and what agrees by construction.  Every number
is read from the records (results/, full_device/full_device_index.json, trace_anchors.json); the text that says what is
tested is fixed here, one entry per figure.
usage: python3 -I evidence_index.py TRACES_DIR [OUT.md]"""
import json
import os
import sys

TR = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(TR, 'EVIDENCE.md')
J = lambda *p: json.load(open(os.path.join(TR, *p)))
A = J('trace_anchors.json')
FD = J('full_device', 'full_device_index.json')
roles = {n: p['role'] for n, p in A['panels'].items()}


def nums(r):
    v1, m, v5 = r['V1_feature_vs_chamfer_px'], r['metric'], r.get('V5_scale_scan') or {}
    nl = r['metric_null']
    worst_null = min(nl['half_pitch_along']['median_abs_um'], nl['5um_across']['median_abs_um'])
    return ('%.2f / %.2f' % (v1['along'], v1['across']), '%.2f / %.2f' % (m['median_abs_um'], m['p90_abs_um']),
            '%.2f' % worst_null,
            ('%.3f / %.3f%s' % (v5['best_along'], v5['best_across'], ' (grid edge)' if v5.get('on_grid_edge') else '')) if v5 else '—')


def fd(name):
    d = FD['devices'].get(name)
    if not d:
        return '—'
    return '`full_device/%s_full.pdf`, `_strip*.png`, `_zoom_*.png` (%.1f × %.1f mm)' % (name, d['footprint_mm'][0], d['footprint_mm'][1])


rows = []


def add(fig, preset, panel_name, rec, tested, consistency):
    v1, med, null, v5 = nums(rec)
    rows.append('| %s | `%s` | `figures/panels/%s.png` | %s | %s | %s | %s | %s | %s | %s | %s |' % (
        fig, preset, panel_name, fd(preset), v1, med, null, v5, roles.get(panel_name, '—'), tested, consistency))


d = J('results', 'dinh', 'dinh_results_rc6.json')
C_D = 'the pitch and the post size (measured on these panels, so V5 along); the scale is the bar\'s'
for k, lab, t, c in (('p1', 'Dinh 2024 Fig. 2B, head', 'the head units\' channel, wall and pitch; the tips', C_D),
                     ('p2', 'Dinh 2024 Fig. 2B, transition', 'the blocks\' length (they read ~2.2 µm shorter: LIMITATIONS 24 a)',
                      C_D + '; the transition\'s form (the blocks beside a 28 µm side) was read on this panel'),
                     ('p3', 'Dinh 2024 Fig. 2B, linear section', 'the side channels\' width at the row the post columns place (%.1f / %.1f '
                      'µm against the schedule\'s %.1f)' % (d['p3']['measured']['image_widths_um']['side_L'], d['p3']['measured']['image_widths_um']['side_R'],
                                                             d['p3']['measured']['mask_widths_um_at_row']['side']),
                      C_D + '; the centre width (the row is the one whose post columns match the image\'s; the outer walls would place it '
                      '%d rows away)' % abs(int(d['p3']['measured']['row_from_post_columns']) - int(d['p3']['measured']['row_from_outer_walls']))),
                     ('p4', 'Dinh 2024 Fig. 2B, outlet', 'the outlet\'s widths against the print (within 0.5 µm; the mask draws the print)',
                      C_D + '; the flare\'s angle and onset (measured on this panel)')):
    add(lab, 'dinh2024_v1_depicted', d[k]['name'], d[k], t, c)
g4 = J('results', 'gifford', 'gifford_results_rc6.json')
add('Gifford 2014 Fig. 4 inset', 'gifford2014_fig4_design', g4['name'], g4,
    'leg 4\'s last rows: centre 100 and sides 300 µm (printed) at the end of the schedule',
    'the pitch and l_g, and the leg spacing (818 µm; read on this inset)')
g6 = J('results', 'gifford', 'gifford6_results_rc6.json')
for j in (7, 5, 3, 1):
    add('Gifford 2014 Fig. 6, leg %d (own map)' % j, 'gifford2014_fig6_design', g6['leg%d' % j]['name'], g6['leg%d' % j],
        'the leg identity; the leg end (END_CHECK); the side channel at the end (reads 5.4–10.0 µm wider: LIMITATIONS 24 e)',
        'the centre (sets the scale), the pitch, l_g')
a6 = J('results', 'gifford', 'gifford6_results_rc6_onemap.json')['four_legs']
t6 = a6['measured']['turns']
add('Gifford 2014 Fig. 6, four legs and their turns (one map)', 'gifford2014_fig6_design', a6['name'], a6,
    'one map for all four legs; where the turns start along the flow (%.1f µm %s of the image\'s)' % (
        abs(t6['centre_offset_mask_minus_image_um']['u']), 'downstream' if t6['centre_offset_mask_minus_image_um']['u'] < 0 else 'upstream'),
    'the leg axes and the turns\' radii (read on this image: walls a median %.1f µm from the mask)' % t6['nearest_mask_edge_um']['median'])
for f, lab, t in (('xia_results_fig3c_rc6.json', 'Xia 2016 Fig. 3c, leg 1\'s first posts', 'the side channels and outer walls at leg 1\'s first rows (V3 channel_across %.1f / %.1f µm; the centre, printed 300 µm, is a scale input)' % tuple(J('results', 'xia', 'xia_results_fig3c_rc6.json')['V3_by_class'][c_]['mean_signed_um'] for c_ in ('channel_across_C', 'channel_across_L'))),
                  ('xia_results_fig3d_rc6.json', 'Xia 2016 Fig. 3d, legs 5 and 3 below the loops', 'the centre widths there (the taper: leg 3 reads 220 against 233 µm) and the legs\' spacing (from Fig. 3a; the loops\' curvature is a scale input)'),
                  ('xia_results_fig3e_rc6.json', 'Xia 2016 Fig. 3e, legs 5 and 3 half way', 'the taper half way (leg 3 reads narrower: LIMITATIONS 24 c)'),
                  ('xia_results_fig3f_rc6.json', 'Xia 2016 Fig. 3f, leg 7 (its own scale)', 'leg 7\'s side channels a quarter of the way down'),
                  ('xia_results_rc6.json', 'Xia 2016 Fig. 3g, outlet', 'the last side channel (624.6 µm against the 600 the printed F draws)')):
    r = J('results', 'xia', f)
    add(lab, 'xia2016_rs7_depicted', r['name'], r, t, 'the printed widths 300 / 150 µm (they set the scale), the pitch and the posts (measured on these panels)')
fa = J('results', 'xia', 'xia_fig3a_device_rc6.json')
sh = fa['V11_group_shifts']
head = ('| paper figure | shipped preset | zoomed trace | whole device | V1 along / across px | median / p90 µm | metric null '
        '(worst) µm | V5 along / across | role | tested | consistency (by construction) |\n|---|---|---|---|---|---|---|---|---|---|---|')
txt = ['# The evidence, figure by figure (rc6)', '',
       'Generated by `code/evidence_index.py` from the trace records; code `%s`.  Each paper figure the traces cover: the '
       'shipped preset drawn over it, the zoomed trace (`figures/panels/`, each panel at its own scale, with `_clean` and '
       '`_residuals` beside it), the whole device (`full_device/`: every polygon of the DXF with the micrographs placed by the '
       'maps their traces stored; V9 checks the round trip, V10 the single scale), and the registration\'s numbers.  V1: the '
       'feature registration against the independent chamfer registration (px); the metric: each mask edge to the nearest '
       'image outline (µm), against its nulls (the mask shifted half a pitch along or 5 µm across; the worse of the two is '
       'shown); V5: the scale the image prefers along / across the flow (1 = the scale used).  *Role*: a claim passes the '
       'limits in `trace_anchors.json` (regression limits, set with these results in view); a recorded panel is reproduced '
       'and reported.  *Tested* and *consistency*: a value read on the same figure it is traced against agrees with it by '
       'construction (PRESET_PROVENANCE, the convention), so only what the trace does not take from the figure is a test.' %
       (A['code_stamp'].get('sha256') if isinstance(A['code_stamp'], dict) else A['code_stamp']), '', head] + rows + [
       '', '**Xia 2016 Fig. 3a, the whole device** (`figures/xia_fig3a_device_rc6.png`, `full_device/xia2016_rs7_depicted_full.pdf`): '
       'posed at its top (rotation %.3f°, chamfer %.2f px there, %.2f px over the whole device), the drawing\'s far loops lie %+.2f mm and its outlet ring %+.2f mm '
       'from the shipped mask\'s along the legs (V11): at the pitch the panels show, Fig. 3a\'s legs hold more rows than the '
       'schedule places (LIMITATIONS 24 d; the variants in `variants/` close it, `README.md`).  V12, each traced panel against '
       'its box on Fig. 3a: %s.' % (fa['pose_top']['rotation_deg'], fa['pose_top']['chamfer_px'], fa['pose_whole']['chamfer_px'],
                                     sh['loops at the bottom']['shift_along_um'] / 1000.0, sh['outlet ring and ports']['shift_along_um'] / 1000.0,
                                     ', '.join('%s %s' % (k, 'inside' if v['ok'] else 'outside (%.1f of the half-size)' % v['offset_over_half_size'])
                                               for k, v in fa['V12_boxes'].items() if 'ok' in v)),
       '', 'Whole-device figures: %s.%s' % (
           '; '.join('`%s` %.1f × %.1f mm (%d micrograph panel%s)' % (n, v['footprint_mm'][0], v['footprint_mm'][1], len(v['panels']),
                                                                  '' if len(v['panels']) == 1 else 's') for n, v in FD['devices'].items()),
           ('  All four shipped devices at one scale: `full_device/full_device_composite.png` (%.0f µm/px; V10 %s).' % (
               FD['composite']['um_per_px'], 'ok' if FD.get('ok') else 'FAILED')) if FD.get('composite') else ''),
       '', 'What the generator was given for each device, and why: `docs/PRESET_PROVENANCE.md` (printed values as printed, '
       'the decoupled targets D and F among them, as App. B and the glossary read them; unprinted values measured on the '
       'figures where they show them).  What the traces do not show: `docs/LIMITATIONS.md` 24–26.', '']
open(OUT, 'w').write('\n'.join(txt))
print('wrote', OUT, len(rows), 'rows')
