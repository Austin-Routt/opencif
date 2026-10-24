"""opencif.published: the papers' devices as Spec objects, with the published values and the assumptions marked."""
from opencif.physics import schedule as S

# Fcl gives Dinh 2024's closed-form outlet-identity F per variant. V1 and V3 are App. B P-C4.2's frozen values
# (glossary \FclosedVone, \FclosedVthree; tests/test_appendix_trace.py). V2's value is in neither App. B nor the
# glossary. It is carried over from the replica generator's device table and is unverified here, and no preset or
# example uses V2.
Fcl = {'V1': 0.5592, 'V2': 0.6169, 'V3': 0.7189}

# The first-generation factories dinh(), xia(), gifford_bead(), and gifford_plt() were removed in RC6, because they
# carried CIF_v2's reference-process posts (g_p 24 to 30, l_g 50.6 to 60) rather than the papers' measured values, and
# no preset, example, test, or tool used them. The replicas below are the ones that tests/test_presets_match_replicas.py
# checks.


# ================================================================================================================
# meander3: replicas of the published devices as specified (the printed parameters, with this work's physics for what
# is not printed) and as depicted (the head drawn with the unit lengths that the paper's figure shows). Unprinted
# dimensions are measured on the papers' own figures (r4/meas_*.py; EXTRACT.md [MEAS]).
# ================================================================================================================
def gifford_fig4():
    """Return the Gifford 2014 Fig. 4 bead device.

    The printed values are f_gap 3.34e-4, a constant w_c of 100, depth 80, and final w_s 300. G is 20.5 (the record's
    measurement on Fig. 4). g_p is 25.3, which is the inset's row pitch of 45.82 um (read at its 300 um bar with the
    bar's anti-aliased ends) minus G (RC6; RC5 had 25.6). l_g is 59.5, the post columns' center-to-center separation
    minus the printed 100 um center width (RC6; RC5 had 56.8, between outline centers). There is no meander, so the side
    channels start closed behind the first post at Gifford's own sealed start (w_s(0) -> 0, p. 4499). The layout
    follows Fig. 4, with a 4-leg paperclip, the inlet inside the second loop, and three ports at the top (two large
    side ports flanking a small center port)."""
    return S.Spec(name='Gifford2014_Fig4', h=80, G=20.5, g_p=R6['gifford']['g_p_um'], l_g=_GIFF_LG, GS=29.7, wc0=100,
                  delta_star=S.x_of(3.34e-4) * 100, F_target=0.90644, ws_in=20.5, W=2400, taper='const',
                  head='none', F_min=4.5, port_r=750.0, port_r_c=300.0, port_r_in=300.0, r_fillet=2.0, rmin_ok=1.8,
                  outlet='separate', pc_last_short=2000.0,
                  note='Gifford 2014 Fig. 4: f 3.34e-4, w_c 100, h 80, final w_s 300; G 20.5 (record), g_p %.1f, '
                       'l_g %.1f (Fig. 4 inset, rc6); sealed start; paperclip 4 legs; separate outlets'
                       % (R6['gifford']['g_p_um'], _GIFF_LG))


def dinh_fig2(variant='V1', f_head=None):
    """Return the Dinh 2024 Fig. 2 device in the given variant.

    The printed values are h 140, GS 28, G 19, w_c 75 -> 52, w_S(in) 21, and f*. g_p is 26.1, the row pitch of
    45.11 um on Fig. 2B's panels 2 to 4 (the 50 um bars read with their anti-aliased ends, 46.91 px) minus G (RC6; RC5
    read 46 px and 27). l_g is 54, from Fig. 2B's post columns, which are 106.06 um apart center to center, minus the
    printed 52. The meander wall is 27 (Fig. 2A inset; Fig. 2B reads 25 to 27). The device has three arms at the chip's
    edge (R0 1281 um, Fig. 2A inset), the depicted inlet and two-stage filter, three rectilinear openings whose side
    channel is GS wide and narrows to w_S(in) at the first post (Fig. 2B, second panel), and one combined filtrate
    outlet around the middle-channel port. With f_head None, the head is drawn as specified (a constant cutoff, with
    three arms and the smallest gap-1 seed that fits W). With f_head given, the head is drawn as depicted (the first
    rows' unit lengths as Fig. 2A's inset draws them)."""
    fstar = {'V1': 1.04e-4, 'V2': 1.12e-4, 'V3': 1.28e-4}[variant]
    tag = 'depicted' if f_head is not None else 'specified'
    return S.Spec(name='Dinh2024_%s_%s' % (variant, tag), h=140, G=19, g_p=R6['dinh']['g_p_um'], l_g=54.0, GS=28.0, wall_M=27.0,
                  wc0=75, delta_star=S.x_of(fstar) * 75, F_target=Fcl[variant], ws_in=21, W=3200.0, margin=15.0, F_min=19.0,
                  taper='shear', wc_min=52, head='dinh3', n_arms_fixed=3, R0_arms=1281.0, stem_half=600.0,
                  filter='dinh', transition='paper', ws_rect=28.0, n_rect=3, outlet='combined',
                  port_r=980.0, port_r_c=_DINH_PORT_C, port_r_in=1020.0, f_head=f_head, seam_at='post',
                  note='Dinh 2024 Fig. 2 (%s head): printed h/GS/G/w_c/w_S(in)/f*; measured g_p %.1f, l_g 54, wall 27, '
                       'R0 1281, filter, ports; F from the frozen outlet identity' % (tag, R6['dinh']['g_p_um']))


# _DINH_PORT_C is the middle-channel port's radius, measured on Dinh's Fig. 2A main panel at 32.3 um/px
# (r4/meas_dinh_layout.py), which replaced an unsourced 830.
import json as _json, os as _os
_DINH_PORT_C = round(_json.load(open(_os.path.join(S.SP, 'r4', 'dinh_layout_meas.json')))['outlet']['centre_port_r_um'], 0)


# R6 holds the figure-derived values as RC6 re-read them, with the scale bars measured with their anti-aliased ends
# and Xia's panels b to g at the scale that their printed widths and loop curvature give (data/r6/rc6_values.json, from
# provenance/rc6_scripts/measure_rc6.py). One rule gives l_g on all three papers, which is the post columns'
# center-to-center separation (free of the outline bias) minus the printed center width at the same place:
#     Dinh     106.06 - 52 = 54.06
#     Gifford  159.5 - 100
#     Xia      220.9 - 150 at the outlet, and 371.0 - 300 at the first posts
# RC5 read l_g between dark outline centers (r4/xia_fig3g_meas.json, 78.1 at the panel-i bar;
# r4/gifford_inset_meas.json, 56.8).
R6 = _json.load(open(_os.path.join(S.SP, 'r6', 'rc6_values.json')))
_XIA_LG = R6['xia']['l_g_um']
_GIFF_LG = R6['gifford']['l_g_um']
_X3A = R6['xia']['layout']


# XIA_TUBE_L is the outlet tubing's length in m, which is assumed, because Xia (p. 10) prints the outlet tubing I.D.s
# (0.25 mm center, 0.86 mm side) without their lengths.
XIA_TUBE_L = 0.15


def tube_R(d, L, mu=1.0e-3):
    """Return the Hagen-Poiseuille resistance of a round tube in SI units, R = 128 mu L / (pi d^4), with mu as in the
    duct law."""
    import math
    return 128.0 * mu * L / (math.pi * d ** 4)


def xia_fig(f_head=None, NM_force=None, n_rect=9):
    """Return the Xia 2016 RS-7 device.

    The printed values are G 16, GS 22.8, w_c(0) 300 -> 150 at constant maximum shear, w_s(1) = G, f_ref 5.8e-4 at
    w_c(ref) 100, depth 135 to 150 (140 used), and about 15:1 filtrate. Panels b to g are read at %.3f um/px, the
    scale that their printed widths (w_c 300 at the first posts and 150 at the outlet), their post columns, and the
    curvature of panel d's loop arcs give (RC6). The single 250 um bar, drawn on panel i, gives 1.936, and RC5 used it
    as 1.953. For the posts (Figs. 3c and 3g), the row pitch is 28.80 px = 51.78 um, so g_p = 51.78 - 16 = 35.8, and
    l_g is 70.9 (the post columns' center-to-center separation minus the printed center width, at the outlet and at
    the first posts). The meander wall is 35.8, from the head's unit pitch (Fig. 3b) and the transition's block pitch
    (Fig. 3c), which are both 65.15 px = 117.1 um = 2 (GS + wall) with the printed GS 22.8. The head envelope W is
    1840 (RC5 had 2000 at the bar). The head has no arms (Figs. 1b(i) and 3b). Nine rectilinear rows follow it (Fig.
    1c(i)), whose side channel is GS wide and narrows to G at the first post (Fig. 1c(ii)). The loops carry Xia's 2 um
    pull-in (p. 4, 'subtract 1-5 um'). RC5 had g_p 40.25, l_g 78.1, and wall 41.

    The layout is taken from Fig. 3a at its 5 mm bar, read with its anti-aliased ends (178.32 px; RC5 had 179;
    r4/meas_xia_fig3a.py x 1.003827). It gives the 7 leg axes, the inlet port (r 862) 2.55 mm downstream of the loop-2
    centers and 1.45 mm above the head's top, the ring outlet (centerline R 4195, with its center 1191 and the
    filtrate port 4781 beyond the center port, and ports r 1723 and 1718), and the ring's top 3.44 mm beyond the
    bottom loops."""
    tag = 'depicted' if f_head is not None else 'specified'
    return S.Spec(name='Xia2016_RS7_%s' % tag, h=140, G=16, g_p=R6['xia']['g_p_um'], l_g=_XIA_LG, GS=22.8,
                  wall_M=R6['xia']['wall_M_um'], wc0=300,
                  delta_star=S.x_of(5.8e-4) * 100, F_target=0.938, ws_in=16, W=R6['xia']['head_width_um'], margin=20.0, F_min=16.0,
                  taper='shear', wc_min=150, head='dinh3', n_arms_fixed=0, stem_half=150.0 + 16.0,
                  transition='paper', ws_rect=22.8, n_rect=n_rect, outlet='combined', outlet_style='ring',
                  port_r=float(_X3A['filtrate_port_radius_um']), port_r_c=float(_X3A['centre_port_radius_um']),
                  port_r_in=float(_X3A['inlet_port_radius_um']), inlet_gap=float(_X3A['inlet_port_gap_um']),
                  ring_R=float(_X3A['ring_radius_um']), ring_dc=float(_X3A['ring_centre_offset_um']),
                  ring_dfp=float(_X3A['ring_filtrate_port_offset_um']), ring_rq=250.0,
                  leg_Y=tuple(float(v) for v in _X3A['leg_Y_um']), pc_inlet_offset=float(_X3A['pc_inlet_offset_um']),
                  pc_out_ext=float(_X3A['pc_out_ext_um']),
                  loop_pull_in=2.0, loop_pull_in_leg=True, loop_ramp=40.0, f_head=f_head, NM_force=NM_force,
                  seam_at='post',
                  R_ext_c=tube_R(0.25e-3, XIA_TUBE_L), R_ext_s=tube_R(0.86e-3, XIA_TUBE_L),
                  note='Xia 2016 RS-7 (%s head): printed G/GS/w_c/w_s(1)/f_ref/F; at %.3f um/px (rc6): posts g_p %.1f, '
                       'l_g %.1f, wall %.1f, head envelope W %.0f; Fig. 3a layout, ports and ring; 2 um loop pull-in'
                       % (tag, R6['xia']['scale_um_per_px'], R6['xia']['g_p_um'], _XIA_LG, R6['xia']['wall_M_um'],
                          R6['xia']['head_width_um']))


if xia_fig.__doc__:                    # the scale, from rc6_values.json (there is no docstring under python -OO)
    xia_fig.__doc__ = xia_fig.__doc__ % R6['xia']['scale_um_per_px']
