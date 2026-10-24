"""Write the full_device.py configuration for rc6's shipped masks and, where built, the two traced variants.
usage: python3 full_device_config.py ROOT RASTERS RESULTS_DIR OUT.json [VARIANTS_DIR]
RESULTS_DIR holds dinh/, gifford/, xia/ (the trace records, tag rc6) and xia/xia_fig3a_device_rc6.json; VARIANTS_DIR
holds the variant builds <name>/<name>.dxf (xia2016_rs7_depicted_dfill, xia2016_rs7_depicted_fig3), whose traces are in
RESULTS_DIR (tags dfill, fig3)."""
import json
import os
import sys

ROOT, RAS, RES, OUT = sys.argv[1:5]
VAR = sys.argv[5] if len(sys.argv) > 5 else None
LAY = json.load(open(os.path.join(ROOT, 'provenance', 'rc6_scripts', 'gifford_layout_measurement.json')))   # the Gifford layouts
run = lambda n, ext: os.path.join(ROOT, 'runs', n, n + ext)
X = os.path.join(RAS, 'xia5-000.png')


def xia_panels(tag):
    out = []
    for L, title in (('c', 'leg 1: the transition\'s last blocks and first posts'),
                     ('d', 'legs 5 and 3 leaving the top loops'), ('e', 'legs 5 and 3 half way down'),
                     ('f', 'leg 7 a quarter of the way down'), ('g', 'the outlet: leg 7\'s last rows and the leads')):
        p = os.path.join(RES, 'xia', ('xia_results_%s.json' % tag) if L == 'g' else ('xia_results_fig3%s_%s.json' % (L, tag)))
        if os.path.exists(p):
            out.append(dict(letter=L, results=p, key='-', raster=X, title='Xia Fig. 3%s, %s' % (L, title)))
    return out


cfg = dict(strip_um_per_px=2.0, composite_um_per_px=5.0, strip_max_mm=40, devices=[
    dict(name='dinh2024_v1_depicted', dxf=run('dinh2024_v1_depicted', '.dxf'),
         title='Dinh 2024 V1 (head as depicted), opencif rc6',
         panels=[dict(letter=k, results=os.path.join(RES, 'dinh', 'dinh_results_rc6.json'), key=k,
                      raster=os.path.join(RAS, 'dinh5-000.png'), title='Dinh Fig. 2B panel %s' % k[1]) for k in ('p1', 'p2', 'p3', 'p4')],
         note=('Every polygon of the shipped DXF; Dinh Fig. 2B\'s four panels at the positions their traces registered '
               '(each panel\'s own 50 µm bar read with its anti-aliased ends).  Dinh prints no whole-device drawing at a '
               'known scale (Fig. 2A\'s main panel read at the 3.2 mm width gives 60.6 mm against the printed 75.8 mm), so '
               'there is no underlay; the device\'s length against the printed 75.8 mm is the conditional conflict of '
               'PRESET_PROVENANCE.')),
    dict(name='gifford2014_fig4_design', dxf=run('gifford2014_fig4_design', '.dxf'),
         title='Gifford 2014 Fig. 4 bead device (design level), opencif rc6',
         panels=[dict(letter='4', results=os.path.join(RES, 'gifford', 'gifford_results_rc6.json'), key='-',
                      raster=os.path.join(RAS, 'gif5-001.png'), title='Gifford Fig. 4 inset, the end of leg 4')],
         note=('Gifford\'s Fig. 4 schematic is not drawn to scale (p. 4502), so there is no underlay; the inset (leg 4\'s '
               'end) is placed where its trace registered.  The preset takes the inset\'s leg spacing (%.0f µm, rc6), so the '
               'neighbouring leg 2 in the inset lies where the image shows it by construction (consistency, not a test): its '
               'posts and centre line on the image\'s, its side walls wider on the image (recorded, not claimed; rc6_traces '
               'README).' % LAY['fig4']['leg_spacing_um'])),
    dict(name='gifford2014_fig6_design', dxf=run('gifford2014_fig6_design', '.dxf'),
         title='Gifford 2014 Fig. 6 platelet device (design level, the micrograph\'s leg axes), opencif rc6',
         panels=[dict(letter='6', results=os.path.join(RES, 'gifford', 'gifford6_results_rc6_onemap.json'), key='four_legs',
                      raster=os.path.join(RAS, 'gif7-000.png'), title='Gifford Fig. 6 micrograph, one map for the four legs and their turns')],
         note=('Gifford\'s Fig. 6 schematic is topology only (not to scale), so there is no underlay; the micrograph (the ends '
               'of legs 7, 5, 3, 1 and the turns of 5, 3, 1) is placed where its one-map trace registered it.  The preset takes '
               'the micrograph\'s layout (rc6): legs 1, 3, 5, 7 on its rows, legs 2, 4, 6 where semicircular turns about the '
               'turn walls\' common centre (%.1f µm beyond leg 1) put them -- only %.0f-%.0f degrees of each turn are in the '
               'frame, so legs 2, 4, 6 are an assumption; leg 7\'s array ends 4 rows beyond the image\'s (LIMITATIONS 25).' % (
                   LAY['fig6']['turn_centre_from_leg1_um'], LAY['fig6']['visible_turn_deg'][0], LAY['fig6']['visible_turn_deg'][1]))),
    dict(name='xia2016_rs7_depicted', dxf=run('xia2016_rs7_depicted', '.dxf'),
         title='Xia 2016 RS-7 (head as depicted), opencif rc6', panels=xia_panels('rc6'),
         underlay=dict(results=os.path.join(RES, 'xia', 'xia_fig3a_device_rc6.json'), key='underlay', raster=X),
         note=('Every polygon of the shipped DXF over Xia\'s Fig. 3a (grey: the paper\'s whole-device drawing, registered '
               'at the device\'s top by xia_fig3a_device.py at its 5 mm bar read with its anti-aliased ends) with Fig. 3\'s '
               'panels c-g at their registered positions (rc6 scale 1.798 µm/px; panel f at its own).  The drawing\'s '
               'loops at the far end and its outlet ring lie 5.7 and 5.4 mm beyond the mask\'s: the schedule places 1079 '
               'rows per full leg against the 1190 Fig. 3a\'s legs hold at the pitch the micrographs show.')),
])
if VAR:
    n = 'xia2016_rs7_depicted_dfill'
    cfg['devices'].append(dict(
        name=n, dxf=os.path.join(VAR, n, n + '.dxf'), title='Xia 2016 RS-7, the D-fill variant (D 6.5315 µm), opencif rc6',
        panels=xia_panels('dfill'),
        underlay=dict(results=os.path.join(RES, 'xia', 'xia_fig3a_device_dfill.json'), key='underlay', raster=X),
        note=('The variant whose rows fill Fig. 3a\'s legs (D 6.5315 µm, F 0.938; provenance/rc6_traces/variants/): the '
              'same drawing and micrographs, registered the same way.')))
    n = 'xia2016_rs7_depicted_fig3'
    if os.path.exists(os.path.join(VAR, n, n + '.dxf')):
        cfg['devices'].append(dict(
            name=n, dxf=os.path.join(VAR, n, n + '.dxf'),
            title='Xia 2016 RS-7 with F and D read off Fig. 3 (F 0.94072, D 6.5705 µm), opencif rc6 variant',
            panels=xia_panels('fig3'),
            underlay=dict(results=os.path.join(RES, 'xia', 'xia_fig3a_device_fig3.json'), key='underlay', raster=X),
            note=('The variant whose last side channel is panel g\'s (624.6 µm: F 0.94072, 15.9 : 1) and whose rows fill '
                  'Fig. 3a\'s legs at that F (D 6.5705 µm; provenance/rc6_traces/code/xia_fig3_fit.py): the same drawing '
                  'and micrographs, registered the same way.')))
json.dump(cfg, open(OUT, 'w'), indent=1)
print('wrote', OUT, [d['name'] for d in cfg['devices']], [len(d['panels']) for d in cfg['devices']])
