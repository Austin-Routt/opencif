# gifford2014_fig6_design_fig6legs

**Verdict: REPRODUCED -- OUTSIDE THE VALIDITY ENVELOPE (u_min 0.0290 below the kappa plateau 0.035): the cutoff law is not calibrated there; cutoff band 3.008-3.624 um, App. B C7 method below the plateau (Hill kappa 4.080 at the widest centre channel) [1-D model: junction band: 49 row(s) beyond +-1%, 12 losing their cutoff (up to 96.1% of D; rows 6032-6080; 0.4% of F)] -- ABOVE THE PROCESS ASPECT-RATIO LIMIT (depth/width 32.0 > 10 at the narrowest straight side channel; reported, not refused or failed in reproduce mode)**

| | target | on the drawn mask |
|---|---|---|
| cutoff D (µm) | 3.6236 | 3.6185 – 3.6287 over all gaps |
| gaps outside D ± 1% (not declared) | 0 | 0 |
| share of the net filtrate extracted within D ± 1% | | 100.000% |
| filtrate fraction F | 0.71768 | 0.717698 |

Design: N = 6080 gaps (0 in the meander head), layout paperclip with 7 leg(s), footprint 51.67 × 5.40 mm.

Cutoff law: D = κ_g δ*, κ_g = 4.9138, 95 % CI [4.750, 5.078] (Gifford 2014 calibration; App. B C7); cutoff band 3.008-3.624 um, App. B C7 method below the plateau (Hill kappa 4.080 at the widest centre channel).

| check | result | value | requirement |
|---|---|---|---|
| DRC at the process floor | ok | `{"lost_open": 0, "lost_solid": 0, "components": 1, "d_um": 4.4}` | nothing narrower than 4.5 um lost (tested at 4.4 um); one fluid component |
| topology | ok | `{"holes": 12159, "expected": 12159}` | one valid polygon; every post and unit present |
| side channels continuous (straight sections, every leg, both sides) | ok | `{"pieces": "all 1", "narrowest_um": 4.692}` | one connected strip per leg and side; narrowest section >= 4.5 um |
| outlet ports separated | ok | `3` | product, filtrate and inlet ports in separate components once the outlet stretch is cut |
| no intersections: parts keep the floor apart | ok | `[{"gap_um": 88.28, "parts": ["inlet and head", "U-turn 2"]}, {"gap_um": 102.982, "parts": ["U-turn 4", "U-turn 6"]}]` | every part (inlet and head, each leg, each U-turn, outlet) at least 4.5 um from every part it is not joined to |
| inlet and outlet ports apart | ok | `46724.9` | every outlet port at least 2000 um (edge to edge) from the inlet port |
| drawn components as requested | ok | `"head none, inlet plain, outlet ring, 7 leg(s)"` | the head style, sealing, inlet filter, outlet style and leg count drawn are the ones the spec asks for (nothing substituted or left empty) |
| head inside its envelope | ok | `62.5` | <= W/2 = 1600 um |
| corner radius | ok | `1.988` | >= 1.80 um |
| every planned feature found on the drawing | ok | `"all"` | every measurement section (head units, gaps, side segments, centre, loops, outlet) finds fluid where the plan puts the feature |
| Q1 drawing realises its design (every gap) | ok | `{"cutoff_dev_max": 0.00142, "at_row": 861, "f_dev_where_no_cutoff": 0.0, "f_row": 1, "F_drawn": 0.717698, "F_intended": 0.717698}` | at every gap the drawn cutoff within 0.25% of the intended one, judged on the larger of the intended cutoff and D (where a row has no cutoff: f within 0.5% of the larger of |f| and f*); F within 1e-4 |
| Q2 filtrate fraction F | ok | `0.717698` | F = 0.71768 +- 0.001 |
| Q2 cutoff D at every gap (undeclared rows) | ok | `{"rows_outside": 0, "first_rows": [], "worst_over": null, "worst_under": null, "undefined_rows": 0, "D_range_um": [3.6185, 3.6287]}` | every gap not covered by a declared option realises D +- 1% on both sides (no reverse flow, f < 1/2) |
| aspect ratio of the narrowest CIF channel (gap openings, centre channel, side channels, U-turn sides, head units, outlet leads, measured; reported, not a pass condition in reproduce mode: LIMITATIONS 7) | **over (reported)** | `{"depth_over_width": 31.97, "at": "straight side channel", "width_um": 4.692}` | depth/width <= process.aspect_ratio_max = 10 (CIF_v2 w_min = max(F_min, h/AR_max) = 15.00 um; spec v5 G7, G11) |
| validity envelope (kappa plateau, gap dominance, Reynolds, Stokes gap and lamina, thin lamina, taper rate) | envelope | `{"u_min": 0.029, "M13": 0.0267, "Re_max": 6.667, "Re_gap_max": 0.00112, "Re_delta_max": 0.00194, "lamina_over_depth": 0.00492, "lamina_over_depth_planned_max": 0.00492, "taper_rate_max": 0.0, "f_design_max": 0.000104, "cutoff_band_um": {"band_um": [3.008380187068954, 3.6236], "basis": "App. B C7 method below the plateau (Hill kappa 4.080 at the widest centre channel)", "retention_side": [3.5028084171109937, 3.744686556229395], "passage_edge_below_knee": 3.008380187068954, "c7_band_constant_lamina": [3.008380187068954, 3.6236], "kappa_at_passage_edge": 4.079528249039471, "method": "below the plateau (u_min < 0.035): App. B P-C7.1, the frozen Hill fit kappa(u) (glossary \\hillKmax, \\hillUhalf, \\hillN) at u = kappa_g delta*/w_c on the widest centre channel; band [kappa(u) delta*, kappa_g delta*], App. B C7's band for a lamina held at delta* (its Xia-scaled reading; for Dinh V1 [\\dinhCollapseLo, \\dinhBandHi] = [1.80, 2.17] um, inside the wider [1.25, 2.17] it records because Dinh's schedule is not printed).  Plateau-anchored, not the implicit fixed point d = kappa(d/w_c) delta, which can lie far lower (App. B C7 CQ3).  retention_side: the kappa_g interval x delta*, the band on the plateau"}}` | u_min >= 0.035 (no gap below the kappa knee); M13 < 0.05; Re <= 10; Re_gap and Re_delta <= 0.01 (CIF_v2 M13/M14); delta(i)/h <= 0.1 at every planned row (E11); |w_c(i) - w_c(i-1)|/(G + g_p) <= 0.01 (M16) |
| 1-D model band (reported, not a pass condition) | reported | `{"junction band (bend, head T-junction and outlet stretch-end terms at both ends of their range)": {"rows_newly_outside_tol": [49, 48], "rows_losing_their_cutoff": 12, "max_cutoff_change_over_D": 0.9608, "rows": "6032-6080", "share_of_F_in_those_rows": 0.004029984143687166}, "t_junction_ignored": {"rows_newly_outside_tol": 0, "rows_losing_their_cutoff": 0, "max_cutoff_change_over_D": 0.0}, "junctions_ignored": {"rows_newly_outside_tol": 58, "rows_losing_their_cutoff": 0, "max_cutoff_change_over_D": 1.9963}, "bend_correction_off": {"rows_newly_outside_tol": 94, "rows_losing_their_cutoff": 0, "max_cutoff_change_over_D": 0.3601}, "outlet balance (both side outlets +-1 %)": {"rows_newly_outside_tol": [48, 48], "rows_losing_their_cutoff": 10, "max_cutoff_change_over_D": 0.9714, "rows": "6033-6080", "share_of_F_in_those_rows": 0.003947329990812895, "dDlast_dlnRside_out": -154.3682}}` | cutoff change on the drawn network when the junction terms (bends, head T-junctions, the outlet stretch's end) take the ends of their range (the depth weight moved; the T-junction term widened; the stretch-end term +-30 %) (and, for information, ignored, and with no bend correction), when both side outlets' resistance moves +-1 % (the outlet balance the last rows rest on), and when the off-chip tubing varies (length x0.5/x2, centre I.D. +-0.025 mm); counted on rows inside D +- tol on the drawn mask |

Status: hard checks PASS; validity envelope OUTSIDE; chip area NOT GIVEN; aspect ratio of CIF channels ABOVE THE LIMIT (reported); length budget G9 NOT GIVEN (no chip); drawn mask assessment RUN: opencif's own checks on the drawn polygons (DRC at the floor, topology, continuity, ports, clearances, corner radii, aspect ratio) and the 1-D network measured on them (Q1, Q2); physics auditor P8 NOT RUN ON THE MASK: CIF_v2's physics auditor (P1-P10) is not run on the exported mask; its P8 criteria are computed from the design as diagnostics, not gated (validation.gates.P8); fabricated or measured NO: Q1 and Q2 are statements of the 1-D network on the drawn mask; reproduction DESIGN-LEVEL: D and F from printed values, drawn with the paper's parts and printed dimensions where printed; not a full-mask match, not a measured-performance match (docs/REPRODUCTIONS.md).

Notes:
- straight side channel w_s = 4.69 um (row 1) below w_min = max(F_min, h/AR_max) = max(4.50, 150.0/10) = 15.00 um (CIF_v2 Process.w_min and design(); spec v5 sec. 18; process.aspect_ratio_max): spec v5 G7 asks widths >= w_min, and G11 would draw that row as a fold, which opencif does not draw: reported, not refused (gates.mode = reproduce)
- measured on the drawing: every side segment, head unit length and width, U-turn side and centre arc, outlet lead, gap opening and centre width.  Taken from the design (not measured): the post length across the flow (gap length l_g) and so the side channels' inner edge, the radii of the head and outlet arcs (their widths are measured), the head rows' centre-path lengths, and any off-chip tubing.
- outside the validity envelope: kappa plateau (u_min 0.0290 < 0.035) (built because gates.mode = reproduce)

Values not given in the spec (defaults and derivations):
- process.gap_entrance_c = 1.0 (default)
- process.aspect_ratio_max = 10.0 (default)
- chip.margin_um = 500.0 (default)
- meander.style = 'pack' (default)
- meander.seal = 'require' (default)
- meander.transition = 'exact' (default)
- meander.rectilinear_rows = 2 (default)
- meander.seam_at = 'segment' (default)
- physics.fail_safe_margin = False (default)
- physics.trim_F = False (default)
- start.style = 'sealed' (default)
- physics.wall_rule = 'exact' (default)
- physics.bend_model = 'stokes' (default)
- meander.channel_um = 1.45 x gap (derived)
- meander.wall_um = max(min_feature, channel/2) (derived)
- meander.seam_side_width_um = gap (derived)
- meander.narrow_um = min_feature (derived)
- meander.head_width_um = the narrowest that seals gap 1, searched from 0.3 to 1 x (chip width - 2 margin) (or 3200) (derived)
- inlet.stem_half_width_um = w_c/2 + max(gap, min_feature) (derived)
