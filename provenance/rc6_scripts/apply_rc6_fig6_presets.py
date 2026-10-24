"""rc6: the two Gifford Fig. 6 presets at the post geometry measured on Fig. 6's own micrograph
(provenance/rc6_scripts/measure_gifford_fig6.py -> gifford_fig6_measurement.json): g_p = the row pitch - G (G 20.5 as
Fig. 4) and l_g; rc5 used Fig. 4's l_g 56.8 and an unscripted pitch reading of 52.7 um.  The values are rounded to
0.1 um from the measurement file (g_p from the unrounded pitch).  Each edit by exact string
replacement (it stops if a target is not found exactly once).  Usage: python3 apply_rc6_fig6_presets.py OPENCIF_ROOT"""
import json
import os
import sys

R = sys.argv[1]
G6 = json.load(open(os.path.join(R, 'provenance', 'rc6_scripts', 'gifford_fig6_measurement.json')))


def edit(rel, old, new, count=1):
    p = os.path.join(R, rel)
    s = open(p).read()
    assert s.count(old) == count, (rel, old[:70], s.count(old))
    open(p, 'w').write(s.replace(old, new))


def fmt(v):
    return ('%.1f' % v).rstrip('0').rstrip('.') if abs(v - round(v)) > 1e-9 else '%d' % round(v)


g6p, g6l = fmt(round(G6['row_pitch_um'] - 20.5, 1)), fmt(round(G6['l_g_um'], 1))
for rel in ('presets/gifford2014_fig6_design.yaml', 'presets/gifford2014_fig6_as_fabricated_15um.yaml'):
    edit(rel, 'post_length_um: 32.2, post_width_um: 56.8,', 'post_length_um: %s, post_width_um: %s,' % (g6p, g6l))
    edit(rel, "  assumed: [G 20.5 and l_g 56.8 as Fig. 4; g_p 32.2 = the inset pitch 52.7 um (measured) - G; ring and port sizes]",
         "  assumed: [G 20.5 as Fig. 4 (the micrograph reads %.1f um along the flow -- relief contrast shifts outlines);\n"
         "    ring and port sizes]\n"
         "  measured_on_figure: [the row pitch %.2f um (%.2f-%.2f; at most %.2f by the caption's 250 um bar which runs\n"
         "    off the raster) and l_g %.1f um (%.1f-%.1f) on Fig. 6's micrograph at its printed 125 um centre --\n"
         "    provenance/rc6_scripts/measure_gifford_fig6.py; g_p = pitch - G; rc5 used Fig. 4's l_g 56.8 and an\n"
         "    unscripted pitch of 52.7 um]"
         % (G6['G_read_um'], G6['row_pitch_um'], G6['row_pitch_range_um'][0], G6['row_pitch_range_um'][1],
            G6['row_pitch_upper_bound_um'], G6['l_g_um'], G6['l_g_range_um'][0], G6['l_g_range_um'][1]))
print('applied to', R)
