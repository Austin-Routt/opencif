"""rc6: the presets and examples at the rc6 figure readings (opencif/data/r6/rc6_values.json, measure_rc6.py; the
Gifford Fig. 6 presets in apply_rc6_fig6_presets.py, after measure_gifford_fig6.py), each edit by exact string replacement (it stops if a
target is not found exactly once).  The Dinh presets' length notes are rewritten separately (apply_dinh_notes.py).
Re-applied on 2026-10-17 after the first rc6 working tree was lost.  Usage: python3 apply_rc6_presets.py OPENCIF_ROOT"""
import json
import os
import sys

R = sys.argv[1]
V = json.load(open(os.path.join(R, 'opencif', 'data', 'r6', 'rc6_values.json')))
dn, gf, xi, lay = V['dinh'], V['gifford'], V['xia'], V['xia']['layout']


def edit(rel, old, new, count=1):
    p = os.path.join(R, rel)
    s = open(p).read()
    assert s.count(old) == count, (rel, old[:70], s.count(old))
    open(p, 'w').write(s.replace(old, new))


def fmt(v):
    return ('%.1f' % v).rstrip('0').rstrip('.') if abs(v - round(v)) > 1e-9 else '%d' % round(v)


# ---- Dinh V1 (both presets) and the two examples on Dinh's process: g_p; the presets' flare onset and E3
for n in ('dinh2024_v1_depicted', 'dinh2024_v1_specified'):
    rel = 'presets/%s.yaml' % n
    edit(rel, 'process: {depth_um: 140, gap_um: 19, post_length_um: 27, post_width_um: 54, min_feature_um: 19}',
         'process: {depth_um: 140, gap_um: 19, post_length_um: %s, post_width_um: 54, min_feature_um: 19}' % fmt(dn['g_p_um']))
    edit(rel, """# outlet (rc5): Fig. 2B, fourth panel: each side channel's outer wall runs straight 67 um from the last gap's centre,
# then flares at 14.5 deg per side (both walls, from their own straight direction: provenance/rc5_scripts/
# measure_dinh_flare.py); the lead width the flare reaches is the generator's (not visible in the panel)""",
         """# outlet (rc5; rc6 onset): Fig. 2B, fourth panel: each side channel's outer wall runs straight %s um from the last
# gap's centre (rc6: rc5's 67 um reading rescaled to the 50 um bar read with its anti-aliased ends, 46.91 px; rc5 46),
# then flares at 14.5 deg per side (both walls, from their own straight direction: provenance/rc5_scripts/
# measure_dinh_flare.py); the lead width the flare reaches is the generator's (not visible in the panel)""" % fmt(dn['flare_start_um']))
    edit(rel, 'flare_deg: 14.5, flare_start_um: 67}', 'flare_deg: 14.5, flare_start_um: %s}' % fmt(dn['flare_start_um']))
    edit(rel, '  source: Dinh et al. 2024 Lab Chip, Fig. 2 and p. 4; measured F',
         '  source: Dinh et al. 2024 Lab Chip, Fig. 2 (p. 917) and p. 916; measured F')
for rel in ('examples/dinh_targets_sealed_head.yaml', 'examples/platelet_serpentine_chip.yaml'):
    edit(rel, 'process: {depth_um: 140, gap_um: 19, post_length_um: 27, post_width_um: 54, min_feature_um: 19}',
         'process: {depth_um: 140, gap_um: 19, post_length_um: %s, post_width_um: 54, min_feature_um: 19}' % fmt(dn['g_p_um']))

# ---- Gifford Fig. 4 and Fig. 5 presets and the examples on Gifford's process: g_p, l_g
gp, lg = fmt(gf['g_p_um']), fmt(gf['l_g_um'])
for rel in ('presets/gifford2014_fig4_as_fabricated_20um.yaml', 'presets/gifford2014_fig4_as_fabricated_8um.yaml',
            'presets/gifford2014_fig5A_as_fabricated_12p5um.yaml', 'presets/gifford2014_fig5B_as_fabricated_12p5um.yaml',
            'presets/gifford2014_fig5C_as_fabricated_12p5um.yaml', 'presets/gifford2014_fig5D_as_fabricated_12p5um.yaml',
            'examples/gifford_physics_linear.yaml', 'examples/gifford_targets_meander_head_19um.yaml'):
    edit(rel, 'post_length_um: 25.6, post_width_um: 56.8,', 'post_length_um: %s, post_width_um: %s,' % (gp, lg))
edit('presets/gifford2014_fig4_design.yaml', """  post_length_um: 25.6    # g_p along the flow
  post_width_um: 56.8     # l_g across the flow (Fig. 4 inset, outline centres)""",
     """  post_length_um: %s    # g_p along the flow (rc6: the inset's row pitch %.2f um - G; rc5 25.6)
  post_width_um: %s     # l_g across the flow (rc6: the post columns' centre-to-centre separation - the printed 100 um centre;
                          # rc5 56.8 between outline centres)""" % (gp, gf['row_pitch_um'], lg))

# ---- Xia RS-7 (all three presets): the panels' scale; Fig. 3a's layout at its corrected bar
for n in ('xia2016_rs7_specified', 'xia2016_rs7_depicted', 'xia2016_rs7_depicted_pullR'):
    rel = 'presets/%s.yaml' % n
    edit(rel, 'process: {depth_um: 140, gap_um: 16, post_length_um: 40.25, post_width_um: 78.1, min_feature_um: 16, min_radius_um: 4.5}',
         'process: {depth_um: 140, gap_um: 16, post_length_um: %s, post_width_um: %s, min_feature_um: 16, min_radius_um: 4.5}'
         % (fmt(xi['g_p_um']), fmt(xi['l_g_um'])))
    edit(rel, 'channel_um: 22.8, wall_um: 41, head_width_um: 2000, seal: allow,',
         'channel_um: 22.8, wall_um: %s, head_width_um: %d, seal: allow,' % (fmt(xi['wall_M_um']), xi['head_width_um']))
    edit(rel, 'inlet: {style: plain, port_radius_um: 859, stem_half_width_um: 166, port_gap_um: 591}',
         'inlet: {style: plain, port_radius_um: %d, stem_half_width_um: 166, port_gap_um: %d}'
         % (lay['inlet_port_radius_um'], lay['inlet_port_gap_um']))
    edit(rel, 'outlet: {style: ring, port_radius_um: 1712, centre_port_radius_um: 1716,',
         'outlet: {style: ring, port_radius_um: %d, centre_port_radius_um: %d,' % (lay['filtrate_port_radius_um'], lay['centre_port_radius_um']))
    edit(rel, 'ring: {radius_um: 4179, centre_offset_um: 1187, filtrate_port_offset_um: 4762, turn_radius_um: 250}}',
         'ring: {radius_um: %d, centre_offset_um: %d, filtrate_port_offset_um: %d, turn_radius_um: 250}}'
         % (lay['ring_radius_um'], lay['ring_centre_offset_um'], lay['ring_filtrate_port_offset_um']))
    edit(rel, '  leg_Y: [0.0, 1327.0, -1739.0, 2221.0, -2640.0, 3254.0, -3959.0]',
         '  leg_Y: [%s]' % ', '.join('%.1f' % v for v in lay['leg_Y_um']))
    edit(rel, '  pc_inlet_offset: 2543.0\n  pc_out_ext: 3428.0',
         '  pc_inlet_offset: %.1f\n  pc_out_ext: %.1f' % (lay['pc_inlet_offset_um'], lay['pc_out_ext_um']))
    edit(rel, '  source: Xia et al. 2016 Sci. Rep., Figs. 1, 3 and p. 4; tubing length 0.15 m assumed (only the I.D.s are printed)',
         '  source: Xia et al. 2016 Sci. Rep., Figs. 1, 3 and p. 4; tubing length 0.15 m assumed (only the I.D.s are printed); '
         'rc6 -- Fig. 3 panels b-g read at %.3f um/px (their printed widths and loop curvature; the panel-i bar gives 1.936), '
         'Fig. 3a at its bar read with its anti-aliased ends (x %.6f)' % (xi['scale_um_per_px'], xi['fig3a_length_factor']))
# the example on Xia's process (refused on M13 by design)
edit('examples/wbc_xia_process_refused.yaml', """# Xia's own process (140 um deep, 16 um gaps, 78 um posts across the flow) with a 200 um centre channel leaves the""",
     """# Xia's own process (140 um deep, 16 um gaps, %s um posts across the flow, rc6) with a 200 um centre channel leaves the""" % fmt(xi['l_g_um']))
edit('examples/wbc_xia_process_refused.yaml', 'process: {depth_um: 140, gap_um: 16, post_length_um: 40, post_width_um: 78, min_feature_um: 16}',
     'process: {depth_um: 140, gap_um: 16, post_length_um: %s, post_width_um: %s, min_feature_um: 16}' % (fmt(xi['g_p_um']), fmt(xi['l_g_um'])))
edit('examples/wbc_xia_process_refused.yaml', 'meander: {enabled: true, style: pack, channel_um: 22.8, wall_um: 41, seal: require}',
     'meander: {enabled: true, style: pack, channel_um: 22.8, wall_um: %s, seal: require}' % fmt(xi['wall_M_um']))
print('applied to', R)
