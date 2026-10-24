"""rc6: opencif/data/r6/figure_anchors.json (hashed; read by tests/test_figure_anchors.py) from rc5's anchors
(opencif/data/r5/figure_anchors.json, kept), rc6's measurements (rc6_measurements.json, measure_rc6.py) and the Gifford
Fig. 4 trace's widths (provenance/rc6_traces/results/gifford/, measured.leg4_widths_um).  The Gifford Fig. 6
micrograph's measurement is not copied here: tests/test_figure_anchors.py reads it from
provenance/rc6_scripts/gifford_fig6_measurement.json (measure_gifford_fig6.py), as it reads the Dinh length notes'
computations from provenance/rc6_scripts/out/.

Usage: python3 make_figure_anchors_r6.py OPENCIF_ROOT"""
import json
import os
import sys

R = sys.argv[1]
H = os.path.join(R, 'provenance', 'rc6_scripts')
a5 = json.load(open(os.path.join(R, 'opencif', 'data', 'r5', 'figure_anchors.json')))
m6 = json.load(open(os.path.join(H, 'rc6_measurements.json')))
gw = json.load(open(os.path.join(R, 'provenance', 'rc6_traces', 'results', 'gifford',
                                 'gifford_results_design.json')))['measured']['leg4_widths_um']
d6, gf, bars = m6['dinh'], m6['gifford'], m6['bars']
out = dict(note='rc6: local anchors measured on the papers\' own figures (publisher rasters), for '
                'tests/test_figure_anchors.py, at each panel\'s scale bar read with its anti-aliased ends '
                '(provenance/rc6_traces/README.md; rc5 read the bars at their full-intensity length, 46 px for Dinh\'s '
                '50 um and 107 px for Gifford\'s 300 um, 2.0 % and 0.5 % short).  Measurements: '
                'provenance/rc6_scripts/rc6_measurements.json (measure_rc6.py) and provenance/rc6_traces/results/.  '
                'Gifford Fig. 6\'s micrograph: provenance/rc6_scripts/gifford_fig6_measurement.json, read directly by '
                'the test.  Global schematic panels carry no scale (not_to_scale).  The rc5 values are kept under '
                'rc5_values.')
p4 = dict(a5['dinh2024_fig2b_panel4'])
p4.update(scale_um_per_px=bars['dinh_fig2B_panel4']['um_per_px'], scale_bar_px=round(bars['dinh_fig2B_panel4']['px_aa'], 4))
p4['flare_deg_per_side'] = dict(a5['dinh2024_fig2b_panel4']['flare_deg_per_side'])
p4['flare_deg_per_side']['basis'] += '; angles are scale-free, unchanged in rc6'
p4['flare_onset_from_last_gap_centre_um'] = dict(left=round(d6['flare_onset_from_last_gap_centre_um']['left'], 1),
                                                 right=round(d6['flare_onset_from_last_gap_centre_um']['right'], 1),
                                                 tolerance=2.0,
                                                 basis='rc5 measure_dinh_flare.py onset rows (px) converted at the corrected bar')
p4['flare_onset_from_last_post_um'] = {k: round(v, 1) for k, v in d6['flare_onset_from_last_post_um'].items()}
p4['post_pitch_um'] = dict(value=round(d6['row_pitch_um'], 2), tolerance=0.2, readings_um=d6['pitch_readings_um'],
                           sd_um=round(d6['row_pitch_sd_um'], 2),
                           basis='pill centres along two post columns on each of panels 2-4 (six readings); tolerance: '
                                 'the bar read to +-0.3 % plus the readings\' spread')
p4['scale_bar_basis'] = 'the bar read with its anti-aliased ends (bar_length_aa; all four panels %.3f px; rc5 46 px)' % (
    bars['dinh_fig2B_panel4']['px_aa'])
out['dinh2024_fig2b_panel4'] = p4
gi = dict(a5['gifford2014_fig4_inset'])
gi.update(scale_um_per_px=gf['scale_um_per_px'], row_pitch_um=round(gf['row_pitch_um'], 2),
          centre_width_um=round(gw['centre'], 1), final_side_width_um=round(0.5 * (gw['side_inner'] + gw['side_outer']), 1),
          scale_bar_px=round(bars['gifford_fig4_inset']['px_aa'], 4), row_pitch_readings_um=gf['row_pitch_readings_um'],
          l_g_um=round(gf['l_g_um'], 1),
          l_g_basis='l_g by the centre-to-centre separation of each leg\'s two post columns minus the printed 100 um '
                    'centre: free of the outline bias and of the relief contrast\'s along-flow shift',
          final_side_widths_um=dict(inner=round(gw['side_inner'], 1), outer=round(gw['side_outer'], 1)),
          widths_basis='leg 4 widths on the registered trace (provenance/rc6_traces/results/gifford/'
                       'gifford_results_design.json, measured.leg4_widths_um); uncertainty 2 px')
out['gifford2014_fig4_inset'] = gi
nts = dict(a5['not_to_scale'])
nts['gifford2014_fig6'] = ('the same sentence: the schematic (left panel), topology only.  The micrograph (right panel) '
                           'is measured at its printed centre width: provenance/rc6_scripts/gifford_fig6_measurement.json')
out['not_to_scale'] = nts
out['rc5_values'] = dict(
    dinh2024_fig2b_panel4={k: a5['dinh2024_fig2b_panel4'][k] for k in ('scale_um_per_px', 'scale_bar_px',
                                                                         'flare_onset_from_last_gap_centre_um',
                                                                         'flare_onset_from_last_post_um', 'post_pitch_um')},
    gifford2014_fig4_inset={k: a5['gifford2014_fig4_inset'][k] for k in ('scale_um_per_px', 'row_pitch_um',
                                                                          'centre_width_um', 'final_side_width_um',
                                                                          'width_pick_uncertainty_um')})
json.dump(out, open(os.path.join(R, 'opencif', 'data', 'r6', 'figure_anchors.json'), 'w'), indent=1)
print('wrote figure_anchors.json')
