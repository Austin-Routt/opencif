# dinh2024_v1_depicted_gp26

**Verdict: REPRODUCED with declared departures (seed: gap 1 above f*, head drawn as depicted, papers' rectilinear transition) -- OUTSIDE THE VALIDITY ENVELOPE (u_min 0.0290 below the kappa plateau 0.035): the cutoff law is not calibrated there; cutoff band 1.805-2.174 um, App. B C7 method below the plateau (Hill kappa 4.080 at the widest centre channel) [1-D model: junction band: 11 row(s) beyond +-1% (up to 28.4% of D; rows 215, 2693-2702; 0.38% of F)]**

| | target | on the drawn mask |
|---|---|---|
| cutoff D (µm) | 2.1742 | 0.9076 – 28.6095 over all gaps |
| gaps outside D ± 1% (not declared) | 0 | 0 |
| share of the net filtrate extracted within D ± 1% | | 96.521% |
| filtrate fraction F | 0.55920 | 0.559283 |

Design: N = 2702 gaps (220 in the meander head), layout linear with 1 leg(s), footprint 148.09 × 3.16 mm.

Cutoff law: D = κ_g δ*, κ_g = 4.9138, 95 % CI [4.750, 5.078] (Gifford 2014 calibration; App. B C7); cutoff band 1.805-2.174 um, App. B C7 method below the plateau (Hill kappa 4.080 at the widest centre channel).

| check | result | value | requirement |
|---|---|---|---|
| DRC at the process floor | ok | `{"lost_open": 0, "lost_solid": 0, "components": 1, "d_um": 18.9}` | nothing narrower than 19.0 um lost (tested at 18.9 um); one fluid component |
| topology | ok | `{"holes": 5479, "expected": 5479}` | one valid polygon; every post and unit present |
| side channels continuous (straight sections, every leg, both sides) | ok | `{"pieces": "all 1", "narrowest_um": 21.016}` | one connected strip per leg and side; narrowest section >= 19.0 um |
| outlet ports separated | ok | `3` | product, filtrate and inlet ports in separate components once the outlet stretch is cut |
| no intersections: parts keep the floor apart | ok | `[]` | every part (inlet and head, each leg, each U-turn, outlet) at least 19.0 um from every part it is not joined to |
| inlet and outlet ports apart | ok | `140823.2` | every outlet port at least 2000 um (edge to edge) from the inlet port |
| drawn components as requested | ok | `"head arms, inlet dinh_filter, outlet teardrop, 1 leg(s)"` | the head style, sealing, inlet filter, outlet style and leg count drawn are the ones the spec asks for (nothing substituted or left empty) |
| head inside its envelope | ok | `1577.7` | <= W/2 = 1600 um |
| corner radius | ok | `4.9` | >= 4.50 um |
| Q1 drawing realises its design (every gap) | ok | `{"cutoff_dev_max": 0.0016, "at_row": 222, "f_dev_where_no_cutoff": 0.00216, "f_row": 229, "F_drawn": 0.559283, "F_intended": 0.559282}` | at every gap the drawn cutoff within 0.25% of the intended one, judged on the larger of the intended cutoff and D (where a row has no cutoff: f within 0.5% of the larger of |f| and f*); F within 1e-4 |
| Q2 filtrate fraction F | ok | `0.559283` | F = 0.55920 +- 0.001 |
| Q2 cutoff D at every gap (undeclared rows) | ok | `{"rows_outside": 0, "first_rows": [], "worst_over": null, "worst_under": null, "undefined_rows": 0, "D_range_um": [0.9076, 28.6095]}` | every gap not covered by a declared option realises D +- 1% on both sides (no reverse flow, f < 1/2) |
| declared departures (reported, not a pass condition) | reported | `{"by_option": {"seed: gap 1 above f*": {"rows_outside": 1, "first_last": [1, 1], "worst_over": 3.658467910812372, "worst_under": 3.658456388155215, "undefined_rows": 0, "share_of_F_oversize": 0.004925447410562054}, "head drawn as depicted": {"rows_outside": 17, "first_last": [2, 18], "worst_over": 1.9696716312046925, "worst_under": 1.6827792993148285, "undefined_rows": 0, "share_of_F_oversize": 0.02249070928374605}, "papers' rectilinear transition": {"rows_outside": 21, "first_last": [216, 236], "worst_over": 13.158911784204818, "worst_under": 0.41744233591033447, "undefined_rows": 7, "share_of_F_oversize": 0.07604154458785344}}, "reverse_flow_share_of_inlet": 0.03943}` | rows each declared option moves, inside its zone, predicted by the design network before measuring |
| aspect ratio of the narrowest CIF channel (gap openings, centre channel, side channels, U-turn sides, head units, outlet leads, measured; reported, not a pass condition in reproduce mode: LIMITATIONS 7) | ok | `{"depth_over_width": 7.37, "at": "gap opening G", "width_um": 19.0}` | depth/width <= process.aspect_ratio_max = 10 (CIF_v2 w_min = max(F_min, h/AR_max) = 19.00 um; spec v5 G7, G11) |
| validity envelope (kappa plateau, gap dominance, Reynolds, Stokes gap and lamina, thin lamina, taper rate) | envelope | `{"u_min": 0.029, "M13": 0.0184, "Re_max": 4.651, "Re_gap_max": 0.000665, "Re_delta_max": 0.00103, "lamina_over_depth": 0.00316, "lamina_over_depth_planned_max": 0.0116, "taper_rate_max": 0.000278, "f_design_max": 0.000216, "cutoff_band_um": {"band_um": [1.8050219505735243, 2.1741549811760246], "basis": "App. B C7 method below the plateau (Hill kappa 4.080 at the widest centre channel)", "retention_side": [2.101680198743562, 2.2468067472041704], "passage_edge_below_knee": 1.8050219505735243, "c7_band_constant_lamina": [1.8050219505735243, 2.1741549811760246], "kappa_at_passage_edge": 4.079523740267386, "method": "below the plateau (u_min < 0.035): App. B P-C7.1, the frozen Hill fit kappa(u) (glossary \\hillKmax, \\hillUhalf, \\hillN) at u = kappa_g delta*/w_c on the widest centre channel; band [kappa(u) delta*, kappa_g delta*], App. B C7's band for a lamina held at delta* (its Xia-scaled reading; for Dinh V1 [\\dinhCollapseLo, \\dinhBandHi] = [1.80, 2.17] um, inside the wider [1.25, 2.17] it records because Dinh's schedule is not printed).  Plateau-anchored, not the implicit fixed point d = kappa(d/w_c) delta, which can lie far lower (App. B C7 CQ3).  retention_side: the kappa_g interval x delta*, the band on the plateau"}}` | u_min >= 0.035 (no gap below the kappa knee); M13 < 0.05; Re <= 10; Re_gap and Re_delta <= 0.01 (CIF_v2 M13/M14); delta(i)/h <= 0.1 at every planned row (E11); |w_c(i) - w_c(i-1)|/(G + g_p) <= 0.01 (M16) |
| 1-D model band (reported, not a pass condition) | reported | `{"junction band (bend, head T-junction and outlet stretch-end terms at both ends of their range)": {"rows_newly_outside_tol": [10, 11], "rows_losing_their_cutoff": 0, "max_cutoff_change_over_D": 0.2842, "rows": "215, 2693-2702", "share_of_F_in_those_rows": 0.003771963391349293}, "t_junction_ignored": {"rows_newly_outside_tol": 197, "rows_losing_their_cutoff": 0, "max_cutoff_change_over_D": 0.0431}, "junctions_ignored": {"rows_newly_outside_tol": 211, "rows_losing_their_cutoff": 0, "max_cutoff_change_over_D": 0.5976}, "bend_correction_off": {"rows_newly_outside_tol": 60, "rows_losing_their_cutoff": 0, "max_cutoff_change_over_D": 0.0128}, "outlet balance (both side outlets +-1 %)": {"rows_newly_outside_tol": [16, 16], "rows_losing_their_cutoff": 4, "max_cutoff_change_over_D": 1.0363, "rows": "2687-2702", "share_of_F_in_those_rows": 0.005468573661271624, "dDlast_dlnRside_out": -169.483}}` | cutoff change on the drawn network when the junction terms (bends, head T-junctions, the outlet stretch's end) take the ends of their range (the depth weight moved; the T-junction term widened; the stretch-end term +-30 %) (and, for information, ignored, and with no bend correction), when both side outlets' resistance moves +-1 % (the outlet balance the last rows rest on), and when the off-chip tubing varies (length x0.5/x2, centre I.D. +-0.025 mm); counted on rows inside D +- tol on the drawn mask |

Status: hard checks PASS; validity envelope OUTSIDE; chip area NOT GIVEN; aspect ratio of CIF channels WITHIN THE LIMIT; length budget G9 NOT GIVEN (no chip); drawn mask assessment RUN: opencif's own checks on the drawn polygons (DRC at the floor, topology, continuity, ports, clearances, corner radii, aspect ratio) and the 1-D network measured on them (Q1, Q2); physics auditor P8 NOT RUN ON THE MASK: CIF_v2's physics auditor (P1-P10) is not run on the exported mask; its P8 criteria are computed from the design as diagnostics, not gated (validation.gates.P8); fabricated or measured NO: Q1 and Q2 are statements of the 1-D network on the drawn mask; reproduction DESIGN-LEVEL: D and F from printed values, drawn with the paper's parts and printed dimensions where printed; not a full-mask match, not a measured-performance match (docs/REPRODUCTIONS.md).

Notes:
- seal: the arms head does not seal gap 1 within the 3200 um head width: gap 1 extracts 13.24 f* (cutoff 7.954 um against D 2.174 um); allowed by meander.seal: allow (gates.mode reproduce), and gap 1 is a declared row
- measured on the drawing: every side segment, head unit length and width, U-turn side and centre arc, outlet lead, gap opening and centre width.  Taken from the design (not measured): the post length across the flow (gap length l_g) and so the side channels' inner edge, the radii of the head and outlet arcs (their widths are measured), the head rows' centre-path lengths, and any off-chip tubing.
- outside the validity envelope: kappa plateau (u_min 0.0290 < 0.035) (built because gates.mode = reproduce)
- rows planned away from f* (declared departures, e.g. an unsealed gap 1) reach Re_gap 8.20e-03 and Re_delta 1.37e-02, above the Stokes gate 0.01: reported, not gated, since those rows claim no cutoff at D (CIF_v2 M14 gates cutoff rows only)
- a row planned away from f* (a declared departure) has planned f_gap 1.38e-03, above kappa_g's calibration cap 5.76e-04: its reported cutoff is extrapolated

Values not given in the spec (defaults and derivations):
- process.fillet_um = 5.0 (default)
- process.gap_entrance_c = 1.0 (default)
- process.aspect_ratio_max = 10.0 (default)
- chip.margin_um = 500.0 (default)
- physics.fail_safe_margin = False (default)
- physics.trim_F = False (default)
- start.style = 'sealed' (default)
- physics.wall_rule = 'exact' (default)
- physics.bend_model = 'stokes' (default)
- process.min_radius_um = min(4.5, 0.25 x min_feature) (derived)
- meander.narrow_um = min_feature (derived)
