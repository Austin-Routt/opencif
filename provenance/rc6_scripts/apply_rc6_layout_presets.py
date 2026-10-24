"""rc6 (after the author's review of 2026-10-17): the Gifford presets take the layout their figures show.  Gifford prints
no layout dimension (the schematics are "not drawn to scale", p. 4502), so rc6 had left the leg spacing to the
generator's rule (723 um for Fig. 6, 969 um for Fig. 4), which the micrograph traces then showed to be wrong: Fig. 6's
legs lie 588 um apart and its turns are arcs about one centre 1745.7 um beyond leg 1 (the uniform paperclip's turns,
radius d/2, are far tighter), and Fig. 4's inset shows legs 2 and 4 817.9 um apart.  The convention (PRESET_PROVENANCE):
a value the paper prints is used as printed; a value it does not print is measured on its figures where they show it,
and left to the generator only where they do not.  The values come from gifford_layout_measurement.json (read on the
figures without the mask); every edit by exact string replacement (it stops if a target is not found exactly once).
Usage: python3 apply_rc6_layout_presets.py OPENCIF_ROOT"""
import json
import os
import sys

R = sys.argv[1]
M = json.load(open(os.path.join(R, 'provenance', 'rc6_scripts', 'gifford_layout_measurement.json')))


def edit(rel, old, new, count=1):
    p = os.path.join(R, rel)
    s = open(p).read()
    assert s.count(old) == count, (rel, old[:70], s.count(old))
    open(p, 'w').write(s.replace(old, new))


# ---- Fig. 6: the leg axes and the last leg's early end (reproduction.leg_Y, .pc_last_short; reproduce mode)
Y = ', '.join(('%g' % v) for v in M['fig6']['leg_Y_um'])
short = round(M['fig6']['last_leg_end_short_um'])
note6 = ("# Layout from the Fig. 6 micrograph (rc6; Gifford prints none): legs 1, 3, 5, 7 on the post rows' tip lines (588 um\n"
         "# apart), legs 2, 4, 6 where semicircular turns about the 16 turn walls' common centre (%.1f um beyond leg 1) put them\n"
         "# -- only 8-23 deg of each turn is in the frame, so that closure is an ASSUMPTION; the last leg's array ends 5 rows\n"
         "# (%d um) before the others.  provenance/rc6_scripts/gifford_layout_measurement.json; LIMITATIONS 25.\n"
         % (M['fig6']['turn_centre_from_leg1_um'], short))
for rel in ('presets/gifford2014_fig6_design.yaml', 'presets/gifford2014_fig6_as_fabricated_15um.yaml'):
    s = open(os.path.join(R, rel)).read()
    first = s.index('name: ')
    edit(rel, s[:first], s[:first] + note6)
    edit(rel, "gates: {mode: reproduce}\n",
         "gates: {mode: reproduce}\nreproduction: {pc_last_short: %d, leg_Y: [%s]}   # the micrograph's layout (rc6)\n" % (short, Y))

# ---- Fig. 4: the leg spacing of the inset (layout.leg_spacing_um; every mode)
d4 = round(M['fig4']['leg_spacing_um'])
for rel in ('presets/gifford2014_fig4_design.yaml', 'presets/gifford2014_fig4_as_fabricated_20um.yaml',
            'presets/gifford2014_fig4_as_fabricated_8um.yaml'):
    edit(rel, 'layout: {style: paperclip, legs: 4}\n',
         'layout: {style: paperclip, legs: 4, leg_spacing_um: %d}   # Fig. 4 inset: legs 2 and 4 %.1f um apart (rc6; '
         'the generator\'s rule gave 969)\n' % (d4, M['fig4']['leg_spacing_um']))
for rel in ('presets/gifford2014_fig4_as_fabricated_20um.yaml', 'presets/gifford2014_fig4_as_fabricated_8um.yaml'):
    edit(rel, 'port radii (schematic not to scale), leg spacing (generator rule)]',
         'port radii (schematic not to scale)]   # the leg spacing: measured on the Fig. 4 inset (rc6)')
print('Fig. 6 leg_Y [%s], pc_last_short %d; Fig. 4 leg spacing %d um' % (Y, short, d4))
