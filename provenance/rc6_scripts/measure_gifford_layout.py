"""rc6: the Gifford layout quantities read on the papers' figures without any mask, gathered into
gifford_layout_measurement.json, which apply_rc6_layout_presets.py writes into the presets.  Every value here is an
image-side reading made by the trace code (the mask enters none of them):
  Fig. 6  leg axes: legs 1, 3, 5, 7 on the post rows' tip lines and the turn walls' common centre
          (provenance/rc6_traces/measurements/gifford_fig6_turns.json, code/gifford6_turns.py); legs 2, 4, 6 at
          2 Y_c - Y of the leg before (semicircular turns, an assumption); the last leg's early end (the image's last
          pills, results/gifford/gifford6_results_rc6.json, measured.last_leg_end_offset_um);
  Fig. 4  the centre lines of legs 2 and 4 on the inset (results/gifford/gifford_results_rc6.json,
          measured.leg_centreline_separation_um.image).
It reads the rc6 trace records, so `apply_all.sh` run on the rc5 tree needs rc6's provenance/rc6_traces/ in place first
(as its README says); the values are read on the images, so this is an order of steps, not a circle.
usage: python3 measure_gifford_layout.py OPENCIF_ROOT [OUT.json]"""
import json
import os
import sys

R = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(R, 'provenance', 'rc6_scripts', 'gifford_layout_measurement.json')
TR = os.path.join(R, 'provenance', 'rc6_traces')
t = json.load(open(os.path.join(TR, 'measurements', 'gifford_fig6_turns.json')))
g6 = json.load(open(os.path.join(TR, 'results', 'gifford', 'gifford6_results_rc6.json')))
g4 = json.load(open(os.path.join(TR, 'results', 'gifford', 'gifford_results_rc6.json')))
own = t['each_wall_own_centre_v_px']
out = dict(
    note=('rc6 (2026-10-17): the layout quantities the papers print none of, read on their own figures without the mask; the '
          'presets take them (provenance/rc6_scripts/apply_rc6_layout_presets.py), so that the shipped masks reproduce the '
          'figures where the papers leave the layout to the generator.  Written by measure_gifford_layout.py.'),
    fig4=dict(leg_spacing_um=g4['measured']['leg_centreline_separation_um']['image'],
              source='Fig. 4 inset (gif5-001.png at its 300 um bar read with its anti-aliased ends): the centre lines of legs 2 and 4 '
                     '(each between its two post columns), provenance/rc6_traces/code/panels_gifford.py; in a 4-leg opencif paperclip '
                     'legs 2 and 4 are neighbours, one spacing apart'),
    fig6=dict(leg_Y_um=t['leg_Y_um'], turn_centre_from_leg1_um=t['centre_from_leg1_axis_um'], turn_residual_sd_um=t['residual_sd_um'],
              turn_centre_sd_um=dict(common_fit=t['centre_sd_px'][1] * t['um_per_px'],
                                     wall_by_wall=own['sd'] * t['um_per_px'],
                                     note='across the legs: the common fit\'s SD, and the scatter of the centres each wall gives on '
                                          'its own (legs 2, 4, 6 at 2 Y_c - Y carry twice the centre\'s error)'),
              visible_turn_deg=[min(t['visible_turn_deg']), max(t['visible_turn_deg'])],
              last_leg_end_short_um=g6['leg7']['measured']['last_leg_end_offset_um'],
              source=('Fig. 6 micrograph (gif7-000.png at its printed 125 um centre, 4.5248 um/px): legs 1, 3, 5, 7 on the post rows\' '
                      'tip lines (provenance/rc6_scripts/gifford_fig6_measurement.json); legs 2, 4, 6 where semicircular turns about '
                      'the 16 turn walls\' common centre put them (provenance/rc6_traces/code/gifford6_turns.py -> '
                      'measurements/gifford_fig6_turns.json; only 8-23 deg of each turn is in the frame, so the semicircle is an '
                      'assumption); the last leg\'s array ends 5 rows (56.8 px) before the others (panels_gifford6.py)')))
json.dump(out, open(OUT, 'w'), indent=1)
print('Fig. 6 leg_Y %s, last leg %.1f um short, centre SD %.1f um (wall by wall %.1f); Fig. 4 spacing %.2f um' % (
    out['fig6']['leg_Y_um'], out['fig6']['last_leg_end_short_um'], out['fig6']['turn_centre_sd_um']['common_fit'],
    out['fig6']['turn_centre_sd_um']['wall_by_wall'], out['fig4']['leg_spacing_um']))
