"""Three questions a reader of the rc6 traces raised (2026-10-17), each answered by one before / after figure drawn from
the records (nothing pasted):

  fig_rc6_fix_gifford6_turns   Gifford Fig. 6: the turns.  Before: the first form (the micrograph's 588 µm leg
                               spacing on the generator's uniform layout: turns of radius 294 µm about a centre one
                               half-spacing beyond leg 1; variants/history588/).  After: the shipped preset, which takes
                               the measured leg axes (the turn walls' common centre beyond leg 1; legs 2, 4, 6 by
                               semicircular turns) -- consistent with the image's turns by construction except along the
                               flow.  Both on one map fitted on the ends of legs 5, 3 and 1; under each, the turn
                               comparison (panels_gifford6.turn_compare).  The numbers are the records'.
  fig_rc6_fix_xia3g_wall       Xia Fig. 3g: the last side channel.  Before: the depicted preset (F 0.938, the printed
                               ~15 : 1: 600 µm).  After: the variant with F and D read off Fig. 3 (624.6 µm).
  fig_rc6_fix_gifford4_start   Gifford Fig. 4: the first side channels.  The design-level preset (sealed start, 4.86 µm
                               at row 1) against the as-fabricated preset (a 20 µm start), the generator's own zoom
                               windows (first rows; the side channel at row 1, both at one scale).  Not a trace.

usage: python3 -I alignment_fixes.py ROOT RASTERS RESULTS_DIR VARIANTS_DIR FIG_DIR [PARTS, default 123]
VARIANTS_DIR: <name>/<name>.dxf, .json for the variants; history588/gifford2014_fig6_design_fig6legs/ for the first one."""
import json
import os
import pickle
import sys
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import composite

ROOT, RAS, RES, VAR, FIG = sys.argv[1:6]
PARTS = sys.argv[6] if len(sys.argv) > 6 else '123'
sys.path[:0] = [ROOT]   # the records' pickles need the engine
os.makedirs(FIG, exist_ok=True)
run = lambda n, ext='.dxf': os.path.join(ROOT, 'runs', n, n + ext)
var = lambda n, ext='.dxf', sub='': os.path.join(VAR, sub, n, n + ext)


def J(*p):
    return json.load(open(os.path.join(*p)))


def turn_lines(rec):
    t = rec['measured']['turns']
    co = t['centre_offset_mask_minus_image_um']
    out = ['turns (%d walls on the image): wall to the mask\'s nearest edge median %.1f µm, p90 %.1f µm' % (
        t['image']['walls'], t['nearest_mask_edge_um']['median'], t['nearest_mask_edge_um']['p90'])]
    if (co['u'] ** 2 + co['v'] ** 2) ** 0.5 <= 50.0:
        # the flow runs along -u (panels_gifford6), so mask - image u < 0 puts the mask's turns downstream of the image's
        out += ['mask\'s turn centre %.1f µm %s, %+.0f µm across from the image\'s;' % (
                    abs(co['u']), 'downstream' if co['u'] < 0 else 'upstream', co['v']),
                'about the image\'s centre the mask\'s walls at radii %+.1f ± %.1f µm from the image\'s' % (
                    t['radii_about_image_centre']['mean_um'], t['radii_about_image_centre']['sd_um'])]
    else:
        out += ['mask\'s turn centre %.0f µm across from the image\'s: not concentric' % abs(co['v'])]
    return out


def last_side(pkl):
    a = pickle.load(open(pkl, 'rb'))['arrays']
    return float(np.asarray(a['ws'])[int(a['N']) - 2])


def part1():
    # ---- (1) Gifford Fig. 6, the turns
    n6 = 'gifford2014_fig6_design_fig6legs'
    b6 = J(RES, 'gifford', 'gifford6_results_fig6legs588.json')['four_legs']
    a6 = J(RES, 'gifford', 'gifford6_results_rc6_onemap.json')['four_legs']
    M6 = J(ROOT, 'provenance', 'rc6_scripts', 'gifford_layout_measurement.json')['fig6']
    t6, ly = a6['measured']['turns'], M6['leg_Y_um']
    du6 = t6['centre_offset_mask_minus_image_um']['u']
    S1 = dict(
        title='Gifford 2014 Fig. 6, the turns', columns=2, width_in=7.0,
        panels=[dict(results=os.path.join(RES, 'gifford', 'gifford6_results_fig6legs588.json'), key='four_legs',
                     raster=os.path.join(RAS, 'gif7-000.png'), dxf=var(n6, sub='history588'), letter='a',
                     title='before: uniform legs 588 µm apart (turn radius 294 µm)', post_max=[31.9, 71.8], lw=0.35,
                     extra_lines=turn_lines(b6)),
                dict(results=os.path.join(RES, 'gifford', 'gifford6_results_rc6_onemap.json'), key='four_legs',
                     raster=os.path.join(RAS, 'gif7-000.png'), dxf=run('gifford2014_fig6_design'), letter='b',
                     title='after: the preset\'s measured leg axes',
                     post_max=[31.9, 71.8], lw=0.35,
                     extra_lines=turn_lines(a6))],
        footer=('Gifford et al. 2014 Lab Chip 14, 4496-4505, Fig. 6 (adapted; © The Royal Society of Chemistry: permission for reproduction to be obtained before publication, README'
                '), with every edge of (a) the first rc6 variant of gifford2014_fig6_design and (b) the shipped rc6 preset '
                'gifford2014_fig6_design, which now takes the micrograph\'s layout (cyan), each on one map fitted on the ends of '
                'legs 5, 3 and 1 (scale 4.5248 µm/px from the printed 125 µm centre; rows at 49.81°).  Both take the micrograph\'s '
                'leg positions for legs 7, 5, 3, 1 and leg 7\'s early end (pc_last_short %d µm).  '
                '(a) The generator\'s uniform paperclip at the micrograph\'s 588 µm spacing: its turns are semicircles of radius '
                'd/2 about a point one half-spacing beyond leg 1 (the generator\'s own rule for this device: 723 µm, radius '
                '362 µm).  (b) The leg axes read on the micrograph (reproduction.leg_Y): the %d turn walls of legs 5, 3, 1 are '
                'concentric arcs about one centre %.1f µm beyond leg 1 (residual SD %.2f µm; the centre\'s SD %.1f µm in the common '
                'fit, %.0f µm between the walls fitted one by one), and semicircular turns about it put legs 2, 4, 6 at %.0f, %.0f, '
                '%.0f µm (2 Y_c - Y, which doubles the centre\'s error); only %.0f-%.0f° of each turn is in the frame, so the '
                'semicircle (and legs 2, 4, 6) is an assumption.  (b) agrees with the image\'s turns by construction (its leg axes '
                'are read on them: consistency, not a test) except along the flow, which nothing read on the turns sets: the mask\'s '
                'turns start %.1f µm %s of the image\'s.  Under each panel: the nearest-outline metric of all claimed edges and, for '
                'the turns alone, the distance from each wall point followed on the image to the mask\'s nearest edge, and the two '
                'sets\' common-centre fits (panels_gifford6.turn_compare).' % (
                    round(M6['last_leg_end_short_um']), t6['image']['walls'], M6['turn_centre_from_leg1_um'], M6['turn_residual_sd_um'],
                    M6['turn_centre_sd_um']['common_fit'], M6['turn_centre_sd_um']['wall_by_wall'], ly[1], ly[3], ly[5],
                    M6['visible_turn_deg'][0], M6['visible_turn_deg'][1], abs(du6), 'downstream' if du6 < 0 else 'upstream')))
    draw('fig_rc6_fix_gifford6_turns', S1)


def draw(stem, S):
    p = os.path.join(FIG, stem + '.json')
    json.dump(S, open(p, 'w'), indent=1)
    composite.main(p, os.path.join(FIG, stem))


def part2():
    # ---- (2) Xia Fig. 3g, the last side channel
    nx = 'xia2016_rs7_depicted_fig3'
    ws0, ws1 = last_side(run('xia2016_rs7_depicted', '.pkl')), last_side(var(nx, '.pkl'))
    g0, g1 = J(RES, 'xia', 'xia_results_rc6.json'), J(RES, 'xia', 'xia_results_fig3.json')
    fit = J(RES, 'xia', 'xia_fig3_fit.json')['result']


    def g_lines(g, ws):
        v5 = g.get('V5_scale_scan') or {}
        return ['last side channel: image %.1f µm, mask %.1f µm; scale fit across %.3f' % (
            g['measured']['left_side_channel_um'], ws, v5.get('best_across', float('nan')))]


    S2 = dict(
        title='Xia 2016 Fig. 3g, the last side channel', columns=2, width_in=7.0,
        panels=[dict(results=os.path.join(RES, 'xia', 'xia_results_rc6.json'), key='-', raster=os.path.join(RAS, 'xia5-000.png'),
                     dxf=run('xia2016_rs7_depicted'), letter='a', title='before: the preset (F 0.938, the printed ~15 : 1)',
                     post_max=[36.8, 71.9], extra_lines=g_lines(g0, ws0)),
                dict(results=os.path.join(RES, 'xia', 'xia_results_fig3.json'), key='-', raster=os.path.join(RAS, 'xia5-000.png'),
                     dxf=var(nx), letter='b', title='after: F %.5f and D %.4f µm read off Figs. 3g and 3a' % (fit['F'], fit['D_um']),
                     post_max=[36.8, 71.9], extra_lines=g_lines(g1, ws1))],
        footer=('Xia et al. 2016 Sci. Rep. 6:35943, Fig. 3g (CC BY 4.0; adapted) at 1.798 µm/px, with the edges of (a) the shipped '
                'opencif rc6 mask xia2016_rs7_depicted and (b) the rc6 variant xia2016_rs7_depicted_fig3 (cyan).  The side '
                'channel\'s width at the end is set by F (the outlet identity): the printed ~15 : 1 (F 0.938) draws %.1f µm against '
                'the image\'s %.1f µm (outline to outline); F %.5f (%.1f : 1) draws %.1f µm, and D %.4f µm (x %.3f of the printed '
                'f_ref\'s) then fills Fig. 3a\'s legs (code/xia_fig3_fit.py).  A reading of the figures under the model\'s premises, '
                'not of the text.' % (ws0, g0['measured']['left_side_channel_um'], fit['F'], fit['filtrate_to_retentate'], ws1,
                                      fit['D_um'], fit['D_over_text'])))
    draw('fig_rc6_fix_xia3g_wall', S2)


def part3():
    # ---- (3) Gifford Fig. 4, the first side channels (the generator's zoom windows; no micrograph shows them)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from opencif import renders as R


    def load_build(n):
        rep = J(run(n, '.json'))
        pk = pickle.load(open(run(n, '.pkl'), 'rb'))
        arr, fluid = pk['arrays'], pk['fluid']
        sc = SimpleNamespace(x=np.asarray(arr['x']), wc=np.asarray(arr['wc']), ws=np.asarray(arr['ws']), kind=np.asarray(arr['kind']),
                             N=int(arr['N']), NM=int(arr['NM']), spec=arr['spec'],
                             info=dict(turns=list(arr['turns']), R_turn=rep['design'].get('R_turn') or 300.0))
        info = dict(legs=arr['legs'], gx=arr['gx'], loops=arr['loops'], inlet=arr['inlet'], outlet=arr['outlet'], x_top=arr['x_top'])
        wins = dict((lab.split(' (')[0], (lab, w)) for lab, w in R.zoom_windows(sc, info, fluid))
        return rep, fluid, wins, float(sc.ws[sc.NM]), float(rep['design']['D_star_um'])


    builds = [('gifford2014_fig4_design', 'design level: sealed start'), ('gifford2014_fig4_as_fabricated_20um', 'as fabricated: 20 µm start')]
    B = [load_build(n) for n, _ in builds]
    # the side-channel window at one scale for both: the larger window, centred on each build's own
    key_sc = [k for k in B[0][2] if k.startswith('side channel at row')][0]
    key_fr = [k for k in B[0][2] if k.startswith('first straight rows')][0]
    size = max(max(b[2][key_sc][1][2] - b[2][key_sc][1][0], b[2][key_sc][1][3] - b[2][key_sc][1][1]) for b in B)
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0))
    for i, ((n, lab), (rep, fluid, wins, ws1, D)) in enumerate(zip(builds, B)):
        w = wins[key_fr][1]
        R.RD.draw(axes[i, 0], fluid, w, style='micro')
        axes[i, 0].set_title('%s  %s: first rows 1-4' % ('ac'[i], lab), fontsize=8, loc='left')
        x0, y0, x1, y1 = wins[key_sc][1]
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        a = 1.6
        w2 = (cx - size / 2, cy - size / a / 2, cx + size / 2, cy + size / a / 2)
        R.RD.draw(axes[i, 1], fluid, w2, style='micro')
        axes[i, 1].set_title('%s  the side channel at row 1: %.2f µm drawn (D %.2f µm)' % ('bd'[i], ws1, D), fontsize=8, loc='left')
        for ax in axes[i]:
            ax.set_xticks([])
            ax.set_yticks([])
            R._scalebar(ax, 0.25)
    import re
    import textwrap
    m_ = re.search(r'depth/width ([0-9.]+) > ([0-9.]+)', B[0][0]['assessment']['verdict'])
    ar_text = ('depth / width %s > %s' % m_.groups()) if m_ else 'see its assessment'
    fig.text(0.01, 0.005, '\n'.join(textwrap.wrap(
        'opencif rc6, Gifford 2014 Fig. 4 bead device, the generator\'s own drawings (fluid light grey, solid dark grey; no micrograph '
        'shows these rows).  (a, b) gifford2014_fig4_design: Gifford\'s sealed-start convention (w_s(0) = 0, Fig. 2; p. 4498) '
        'starts the side channel from zero width, hence the pointed start; at row 1 it is %.2f µm wide, narrower than the '
        'particles the first gaps are meant to pass (up to D = %.2f µm).  The record reports it above the aspect-ratio limit '
        '(%s); it is a picture of the design equations, not a device to make (LIMITATIONS 26).  (c, d) gifford2014_fig4_as_fabricated_20um: the side channel starts at 20 µm (an '
        'assumed floor; Gifford prints none) and widens from there.  b and d at one scale.' % (B[0][3], B[0][4], ar_text), 175)),
        fontsize=6.5, ha='left', va='bottom')
    fig.tight_layout(rect=(0, 0.09, 1, 1), h_pad=2.2)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(FIG, 'fig_rc6_fix_gifford4_start.' + ext), dpi=300)
    print('wrote the three figures in', FIG)


if '1' in PARTS:
    part1()
if '2' in PARTS:
    part2()
if '3' in PARTS:
    part3()
