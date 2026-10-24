"""rc6: engine/devspecs.py, the replica factories (tests/test_presets_match_replicas.py compares five presets with
them) at the rc6 values read from data/r6/rc6_values.json, and the unused first-generation factories removed.  Each
edit by exact string replacement (it stops if a target is not found exactly once).  Re-applied on 2026-10-17 after the
first rc6 working tree was lost.  Usage: python3 apply_rc6_devspecs.py OPENCIF_ROOT"""
import os
import sys

P = os.path.join(sys.argv[1], 'opencif', 'engine', 'devspecs.py')
s = open(P).read()
a = s.index("def dinh(variant='V1', head='pack'")
b = s.index("# ================================================================================================================\n# meander3:")
s = s[:a] + """# rc6: the first-generation factories dinh(), xia(), gifford_bead() and gifford_plt() are removed: they carried
# CIF_v2's reference-process posts (g_p 24-30, l_g 50.6-60), not the papers' measured values, and no preset,
# example, test or tool used them.  The replicas below are the ones tests/test_presets_match_replicas.py checks.


""" + s[b:]
REP = [
    ("""    G 20.5 / g_p 25.6 (the record's measurement on Fig. 4), l_g 56.8 (Fig. 4 inset, between outline centres;
    r4/gifford_inset_meas.json).  No meander: the side""",
     """    G 20.5 (the record's measurement on Fig. 4) / g_p 25.3 (rc6: the inset's row pitch 45.82 um at its 300 um bar read
    with its anti-aliased ends, minus G; rc5 25.6), l_g 59.5 (rc6: the post columns' centre-to-centre separation minus
    the printed 100 um centre; rc5 56.8 between outline centres).  No meander: the side"""),
    ("""    return S.Spec(name='Gifford2014_Fig4', h=80, G=20.5, g_p=25.6, l_g=_GIFF_LG, GS=29.7, wc0=100,""",
     """    return S.Spec(name='Gifford2014_Fig4', h=80, G=20.5, g_p=R6['gifford']['g_p_um'], l_g=_GIFF_LG, GS=29.7, wc0=100,"""),
    ("""                  note='Gifford 2014 Fig. 4: f 3.34e-4, w_c 100, h 80, final w_s 300; G/g_p 20.5/25.6 (record), '
                       'l_g 56.8 (Fig. 4 inset, outline centres); sealed start; paperclip 4 legs; separate outlets')""",
     """                  note='Gifford 2014 Fig. 4: f 3.34e-4, w_c 100, h 80, final w_s 300; G 20.5 (record), g_p %.1f, '
                       'l_g %.1f (Fig. 4 inset, rc6); sealed start; paperclip 4 legs; separate outlets'
                       % (R6['gifford']['g_p_um'], _GIFF_LG))"""),
    ("""    \"\"\"Dinh 2024 Fig. 2 device: h 140, GS 28, G 19, w_c 75 -> 52, w_S(in) 21, f* (printed); g_p 27, l_g 54 (Fig. 2B
    micrographs, 50 um bars); meander wall 27""",
     """    \"\"\"Dinh 2024 Fig. 2 device: h 140, GS 28, G 19, w_c 75 -> 52, w_S(in) 21, f* (printed); g_p 26.1 (rc6: the row
    pitch 45.11 um on Fig. 2B's panels 2-4, the 50 um bars read with their anti-aliased ends, 46.91 px, minus G; rc5 read
    46 px and 27), l_g 54 (Fig. 2B: the post columns 106.06 um apart centre to centre, minus the printed 52); meander wall 27"""),
    ("""    return S.Spec(name='Dinh2024_%s_%s' % (variant, tag), h=140, G=19, g_p=27.0, l_g=54.0, GS=28.0, wall_M=27.0,""",
     """    return S.Spec(name='Dinh2024_%s_%s' % (variant, tag), h=140, G=19, g_p=R6['dinh']['g_p_um'], l_g=54.0, GS=28.0, wall_M=27.0,"""),
    ("""                  note='Dinh 2024 Fig. 2 (%s head): printed h/GS/G/w_c/w_S(in)/f*; measured g_p 27, l_g 54, wall 27, '
                       'R0 1281, filter, ports; F from the frozen outlet identity' % tag)""",
     """                  note='Dinh 2024 Fig. 2 (%s head): printed h/GS/G/w_c/w_S(in)/f*; measured g_p %.1f, l_g 54, wall 27, '
                       'R0 1281, filter, ports; F from the frozen outlet identity' % (tag, R6['dinh']['g_p_um']))"""),
    ("""# post lengths across the flow, measured between dark outline centres (the convention that reproduces the printed
# channel widths on the same micrographs): Xia Fig. 3g (r4/meas_xia_fig3g.py; was 84, an outer-edge reading) and
# Gifford's Fig. 4 inset (r4/meas_gifford_inset.py; was 60)
_XIA_LG = round(float(np.mean(_json.load(open(_os.path.join(S.SP, 'r4', 'xia_fig3g_meas.json')))['post_across_flow_um'])), 1)
_GIFF_LG = round(float(np.mean(_json.load(open(_os.path.join(S.SP, 'r4', 'gifford_inset_meas.json')))['post_across_flow_um_line_centres'])), 1)""",
     """# rc6: the figure-derived values re-read with the scale bars measured with their anti-aliased ends, and Xia's panels
# b-g at the scale their printed widths and loop curvature give (data/r6/rc6_values.json, from
# provenance/rc6_scripts/measure_rc6.py).  l_g by one rule on all three papers: the post columns' centre-to-centre
# separation (free of the outline bias) minus the printed centre width at the same place (Dinh 106.06 - 52 = 54.06;
# Gifford 159.5 - 100; Xia 220.9 - 150 at the outlet, 371.0 - 300 at the first posts).  rc5 read l_g between dark
# outline centres (r4/xia_fig3g_meas.json 78.1 at the panel-i bar; r4/gifford_inset_meas.json 56.8).
R6 = _json.load(open(_os.path.join(S.SP, 'r6', 'rc6_values.json')))
_XIA_LG = R6['xia']['l_g_um']
_GIFF_LG = R6['gifford']['l_g_um']
_X3A = R6['xia']['layout']"""),
    ("""    Posts from the Fig. 3g micrograph (250 um bar; w_c reads 150 um): pitch 56.25 um (10-period autocorrelation of
    both post rows) -> g_p = 56.25 - 16 = 40.25; l_g 78.1 (across-flow, outline centres, r4/xia_fig3g_meas.json).  (Fig. 1c(ii) on its label scale gave
    37.4 / 76; Fig. 1 is a schematic whose bars disagree with its labels, so the micrograph is used.)
    Meander wall 41 and head envelope from the Fig. 3b micrograph (row pitch 127.5 um = 2 (GS + wall), first tip
    939 um from the axis); no arms (Fig. 1b(i), 3b); nine rectilinear rows (Fig. 1c(i)) whose side channel is GS wide
    and narrows to G at the first post (Fig. 1c(ii)); Xia's 2 um loop pull-in (p. 4: 'subtract 1-5 um').
    Layout from Fig. 3a at its 5 mm bar (179 px, 27.93 um/px; r4/meas_xia_fig3a.py): the 7 leg axes, the inlet port
    (r 859) 2.54 mm downstream of the loop-2 centres and 1.45 mm above the head's top, the ring outlet (centreline R 4179, centre 1187 and filtrate
    port 4762 beyond the centre port; ports r 1716 / 1712), the ring's top 3.43 mm beyond the bottom loops.\"\"\"""",
     """    rc6: panels b-g are read at %.3f um/px, the scale their printed widths (w_c 300 at the first posts, 150 at the
    outlet), their post columns and the curvature of panel d's loop arcs give (the single 250 um bar, drawn on panel i,
    gives 1.936; rc5 used it as 1.953).  Posts (Figs. 3c, 3g): row pitch 28.80 px = 51.78 um -> g_p = 51.78 - 16 = 35.8;
    l_g 70.9 (the post columns' centre-to-centre separation minus the printed centre width, at the outlet and at the
    first posts).  Meander wall 35.8 from the head's unit pitch (Fig. 3b) and the transition's block pitch (Fig. 3c),
    both 65.15 px = 117.1 um = 2 (GS + wall) with the printed GS 22.8; head envelope W 1840 (rc5's 2000 at the bar);
    no arms (Fig. 1b(i), 3b); nine rectilinear rows (Fig. 1c(i)) whose side channel is GS wide and narrows to G at the
    first post (Fig. 1c(ii)); Xia's 2 um loop pull-in (p. 4: 'subtract 1-5 um').  (rc5: g_p 40.25, l_g 78.1, wall 41.)
    Layout from Fig. 3a at its 5 mm bar read with its anti-aliased ends (178.32 px; rc5 179; r4/meas_xia_fig3a.py x
    1.003827): the 7 leg axes, the inlet port (r 862) 2.55 mm downstream of the loop-2 centres and 1.45 mm above the
    head's top, the ring outlet (centreline R 4195, centre 1191 and filtrate port 4781 beyond the centre port; ports r
    1723 / 1718), the ring's top 3.44 mm beyond the bottom loops.\"\"\""""),
    ("""    return S.Spec(name='Xia2016_RS7_%s' % tag, h=140, G=16, g_p=40.25, l_g=_XIA_LG, GS=22.8, wall_M=41.0, wc0=300,
                  delta_star=S.x_of(5.8e-4) * 100, F_target=0.938, ws_in=16, W=2000.0, margin=20.0, F_min=16.0,""",
     """    return S.Spec(name='Xia2016_RS7_%s' % tag, h=140, G=16, g_p=R6['xia']['g_p_um'], l_g=_XIA_LG, GS=22.8,
                  wall_M=R6['xia']['wall_M_um'], wc0=300,
                  delta_star=S.x_of(5.8e-4) * 100, F_target=0.938, ws_in=16, W=R6['xia']['head_width_um'], margin=20.0, F_min=16.0,"""),
    ("""                  port_r=1712.0, port_r_c=1716.0, port_r_in=859.0, inlet_gap=1450.0 - 859.0,
                  ring_R=4179.0, ring_dc=1187.0, ring_dfp=4762.0, ring_rq=250.0,
                  leg_Y=(0.0, 1327.0, -1739.0, 2221.0, -2640.0, 3254.0, -3959.0), pc_inlet_offset=2543.0,
                  pc_out_ext=3428.0,""",
     """                  port_r=float(_X3A['filtrate_port_radius_um']), port_r_c=float(_X3A['centre_port_radius_um']),
                  port_r_in=float(_X3A['inlet_port_radius_um']), inlet_gap=float(_X3A['inlet_port_gap_um']),
                  ring_R=float(_X3A['ring_radius_um']), ring_dc=float(_X3A['ring_centre_offset_um']),
                  ring_dfp=float(_X3A['ring_filtrate_port_offset_um']), ring_rq=250.0,
                  leg_Y=tuple(float(v) for v in _X3A['leg_Y_um']), pc_inlet_offset=float(_X3A['pc_inlet_offset_um']),
                  pc_out_ext=float(_X3A['pc_out_ext_um']),"""),
    ("""                  note='Xia 2016 RS-7 (%s head): printed G/GS/w_c/w_s(1)/f_ref/F; Fig. 3g posts g_p 40.25, l_g 78.1; '
                       'Fig. 3b wall 41 and head envelope W 2.0 mm; Fig. 3a layout, ports and ring; 2 um loop pull-in'
                       % tag)""",
     """                  note='Xia 2016 RS-7 (%s head): printed G/GS/w_c/w_s(1)/f_ref/F; at %.3f um/px (rc6): posts g_p %.1f, '
                       'l_g %.1f, wall %.1f, head envelope W %.0f; Fig. 3a layout, ports and ring; 2 um loop pull-in'
                       % (tag, R6['xia']['scale_um_per_px'], R6['xia']['g_p_um'], _XIA_LG, R6['xia']['wall_M_um'],
                          R6['xia']['head_width_um']))


xia_fig.__doc__ = xia_fig.__doc__ % R6['xia']['scale_um_per_px']     # rc6: the scale, from rc6_values.json"""),
]
for old, new in REP:
    assert s.count(old) == 1, (old[:80], s.count(old))
    s = s.replace(old, new)
open(P, 'w').write(s)
print('applied to', P)
