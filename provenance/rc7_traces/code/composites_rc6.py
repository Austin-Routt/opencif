"""The rc6 dissertation figures redrawn from the trace records (composite.py): Dinh Fig. 2B (the four panels), Gifford's
Fig. 4 inset, Xia's Fig. 3 panels c-g (each at its registered map; Gifford's Fig. 6 figures are panels_gifford6.py's
own whole-panel figures).  usage: python3 composites_rc6.py ROOT RASTERS RESULTS_DIR FIG_DIR"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import composite

ROOT, RAS, RES, FIG = sys.argv[1:5]
os.makedirs(FIG, exist_ok=True)
run = lambda n: os.path.join(ROOT, 'runs', n, n + '.dxf')
G4 = json.load(open(os.path.join(RES, 'gifford', 'gifford_results_rc6.json')))
# the composite colours leg 4 (the claim) by the mid line between the legs; the trace classifies by each edge's fluid side
# and records where the two differ -- the colouring is valid only if they do not (panels_gifford.py, claim_split)
assert G4['measured'].get('claim_split', {}).get('disagree_with_mid_line') == 0, 'Fig. 4: the claim split differs from the mid line'
specs = {
    'fig_dinh2024_fig2B_traced_rc6': dict(
        title='Dinh 2024 Fig. 2B', columns=4, width_in=7.0,
        panels=[dict(results=os.path.join(RES, 'dinh', 'dinh_results_rc6.json'), key=k, raster=os.path.join(RAS, 'dinh5-000.png'),
                     dxf=run('dinh2024_v1_depicted'), letter=L, title=t, post_max=[27.1, 55.0])
                for k, L, t in (('p1', 'a', 'head'), ('p2', 'b', 'transition'), ('p3', 'c', 'linear section'), ('p4', 'd', 'outlet'))],
        footer=('Dinh et al. 2024 Lab Chip 24, 913-923, Fig. 2B (adapted; © The Royal Society of Chemistry: permission for reproduction to be obtained before publication, README'
                '), with the edges of the shipped opencif rc6 mask dinh2024_v1_depicted (cyan).  Scale: each panel\'s own '
                '50 µm bar read with its anti-aliased ends (46.91 px; 1.0660 µm/px); rotation measured on the image (0.36-0.63°); '
                'translation fitted on features and, independently, by chamfer matching (they agree within 1.1 px).  Under '
                'each panel: the median and 90th percentile of the distance from the mask\'s edges to the image\'s outlines, the '
                'same with the mask shifted half a row pitch along and 5 µm across (the nulls), and the scale the image prefers '
                'along and across the flow (1 = the bar\'s).  Row pitch Λ = G + g_p = 19 + 26.1 = 45.1 µm; link length ℓ_g 54 µm.')),
    'fig_gifford2014_fig4_inset_traced_rc6': dict(
        title='Gifford 2014 Fig. 4 inset', columns=1, width_in=3.4,
        panels=[dict(results=os.path.join(RES, 'gifford', 'gifford_results_rc6.json'), key='-', raster=os.path.join(RAS, 'gif5-001.png'),
                     dxf=run('gifford2014_fig4_design'), letter='',
                     title=('the end of leg 4 (cyan) beside the neighbouring leg at the inset\'s spacing (amber)'
                            if G4['measured'].get('spacing_from_this_inset') else 'the end of leg 4 (cyan) beside the neighbouring leg (amber)'),
                     claim_y_above=0.5 * (G4['mask_record']['leg_Y'][1] + G4['mask_record']['leg_Y'][3]),
                     post_max=[26.3, 60.5])],
        footer=('Gifford et al. 2014 Lab Chip 14, 4496-4505, Fig. 4 inset (adapted; © The Royal Society of Chemistry: permission for reproduction to be obtained before publication, README'
                '), with the shipped opencif rc6 mask gifford2014_fig4_design: leg 4\'s last rows (cyan, the '
                'claim) and the neighbouring leg at the preset\'s leg spacing (amber; Gifford prints none, rc6 reads it on this inset, '
                '818 µm, so its centre line is the image\'s by construction; its side channels\' widths there depend on which of its '
                'rows pass the inset).  '
                'Scale: the 300 µm bar read with its anti-aliased ends (107.52 px); the legs at 45.78° on the image.')),
    'fig_xia2016_fig3_traced_rc6': dict(
        title='Xia 2016 Fig. 3 c-g', columns=5, width_in=7.0,
        panels=[dict(results=os.path.join(RES, 'xia', f), key='-', raster=os.path.join(RAS, 'xia5-000.png'),
                     dxf=run('xia2016_rs7_depicted'), letter=L, title=t, post_max=[36.8, 71.9])
                for f, L, t in (('xia_results_fig3c_rc6.json', 'c', 'leg 1, first posts'), ('xia_results_fig3d_rc6.json', 'd', 'legs 5, 3 below the loops'),
                                ('xia_results_fig3e_rc6.json', 'e', 'legs 5, 3 half way'), ('xia_results_fig3f_rc6.json', 'f', 'leg 7'),
                                ('xia_results_rc6.json', 'g', 'the outlet'))],
        footer=('Xia et al. 2016 Sci. Rep. 6:35943, Fig. 3 c-g (CC BY 4.0; adapted), with the shipped opencif rc6 mask '
                'xia2016_rs7_depicted (cyan).  Scale: 1.798 µm/px, the scale the panels\' printed widths, post columns and loop '
                'curvature give (the only bar is panel i\'s); panel f at its own pitch (1.774 µm/px).  Panel e is recorded, not '
                'claimed: leg 3 there is narrower than the constant-shear taper puts it (186 vs 205 µm), so one translation '
                'cannot fit both legs.')),
}
for stem, S in specs.items():
    p = os.path.join(FIG, stem + '.json')
    json.dump(S, open(p, 'w'), indent=1)
    composite.main(p, os.path.join(FIG, stem))
