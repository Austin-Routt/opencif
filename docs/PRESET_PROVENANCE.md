# Where every preset value comes from

Every preset value carries one of five codes. P means printed in the paper. M means measured on the paper's figure,
which is an image measurement, a reading of the published raster at its bar or at a printed width rather than a
measurement of a device (the scripts are in the replica package `r4/`). In rc6, every figure-derived value was re-read
with each scale bar read with its anti-aliased ends and with two uniform rules, l_g = the post columns'
center-to-center separation minus the printed center, and g_p = the row pitch minus G (`provenance/rc6_scripts/`, with
the values in `opencif/data/r6/rc6_values.json`). D means derived from printed values by the frozen framework. Since
rc7, D (fitted) means derived through the framework from a figure measurement, so such a value reproduces that
measurement by construction and is tested only by the others. A means assumed, because the paper does not give the
value. G means chosen by the generator or by us to make the build work, which makes the value a fit rather than a fact
about the device. The reader should treat the values marked A or G as open. Rows marked *departure* deliberately draw
a row off the constant cutoff and are allowed only in `gates.mode: reproduce`.

The convention (rc6, after the author's review of 2026-10-17) is that a value the paper prints is used as printed (the
decoupled targets D and F among them, as App. B and the glossary read them), whereas a value that it does not print is
measured on its figures where they show it (code M) and is left to the generator (code G) or assumed (code A) only
where they do not. The layout is the main case. None of the three papers prints one, so the Gifford presets take the
leg axes and spacing that their micrographs show (`provenance/rc6_scripts/gifford_layout_measurement.json`,
`apply_rc6_layout_presets.py`), as Xia's presets take Fig. 3a's. Where a figure contradicts a printed value (Xia's last
side channel on Fig. 3g against the printed ~15 : 1), the preset keeps the printed value and the figure's reading is a
traced variant (`provenance/rc6_traces/variants/`), so the difference is reported rather than absorbed. Two values that
a figure shows stay A, each with its reason in its row. First, Gifford Fig. 6's gap G keeps Fig. 4's 20.5 µm, because
the micrograph's 17.1 µm is a length along the flow between two outlines, which relief contrast shifts by a few µm each
(the Fig. 4 inset reads its outlines long along the flow), so the reading is reported and is not trusted to set G.
Second, Xia's center-channel profile between its printed ends is opencif's constant-shear form of Xia's printed
principle (a relatively constant maximum shear along each leg, found numerically; pp. 3–4), because the generator takes
no user taper profile (LIMITATIONS 16), so the panels' widths along the taper are reported (LIMITATIONS 24 c) rather
than drawn. The aligned presets of rc7 do what this paragraph's last two cases do not. They draw the gaps that the
micrographs read (Gifford Figs. 4 and 6, G 17.2 and 17.1 µm, code M, at the rc7 trace window) and the taper that Xia's
panels read (the parallel-plate form, a reproduce-mode departure), whereas the design-level presets keep the rule
above. Values read on a figure reproduce that figure by construction, so their agreement with it is consistency rather
than an independent test.

## Framework constants (all presets)

| value | code | source / note |
|---|---|---|
| cutoff law D = κ_g x(f) w_c, cubic kernel f = 3x² − 2x³ | D | frozen v114 record; κ_g = 4.9138 [4.750, 5.078] (±3.4%) fitted on Gifford 2014 Fig. 3 |
| below the plateau: κ(u) Hill fit, κ_max 5.568, u_½ 0.01651, n 1.791 | D | glossary macros; App. B P-C7.1's method gives the reported cutoff band |
| validity plateau u = D/w_c ≥ 0.035, foot 0.012 | A | external inputs (glossary, Validity Envelope; App. B C8 CQ7, "its value is not derived"), which are an empirical bin edge in the empty interval 0.0294–0.0353 of Gifford Fig. 3's points and are not tuned |
| gap-dominance gate M13 < 0.05 (spec v5 §19, CIF_v2 M13); channel Re ≤ 10 and Stokes gates Re_gap, Re_δ ≤ 10⁻² (CIF_v2 M14; App. B thesis scope), water | D | gates read from `CIF_v2.Tolerances` |
| duct law r(w, h) (Gifford eqn 5, 50 terms) | P | Gifford 2014 eqn 5 |
| gap resistance r(G)(l_g + c_ent D_h), c_ent = 1 | A | spec §14.4 form; the entrance coefficient is not measured |
| curved-duct correction on every arc (U-turns, head hairpins and arms, transition unit, outlet leads) | D | the exact fully developed Stokes solution for a rectangular duct curved in its plane (Fourier-Bessel series, `schedule.bend_ratio`), tested against its planar-annulus and Hele-Shaw limits and an independent finite-difference solve; Dean inertia (≤ 1.5e-4 at Re ≤ 10, a/R ≤ 0.2) is left out |
| bend junction coefficient κ_j: 0.50 (deep) … 0.33 (shallow) | G | a fit to 2-D Stokes solves of straight-arc-straight channels in the deep and shallow limits at four geometries (centerline radius / width 0.98–4.8; reviews R3 and R4, `provenance/junction_fits/r3_junction2d*.py`; deep values 0.46–0.52); each arc's two junctions remove κ_j·w from the developed correction |
| head T-junction c_T: −0.28 w (deep) … −0.47 w (shallow) | G | a fit to 2-D solves of a unit's drop turning into the next unit's stem with the gap's mouth below (review R4, `provenance/junction_fits/r4_tjunc2d*.py`; plain 90° corner −0.24 / −0.44, the latter the textbook Hele-Shaw corner), added to every taper and arm unit |
| head T-junction fillet term (r_f/w) × S_T: S_T −0.83 (deep) … −0.33 (shallow), r_f/w capped at 0.3 | G | a fit to 2-D solves of the filleted T at r_f/w 0.18 and 0.26 (review R4b, `provenance/junction_fits/r4b_fillet*.py`, with deep −0.143 / −0.21 ± 0.03 and shallow −0.055 / −0.088); at the heads' widths it is about a third of the T term (28 µm units, 5 µm fillets: −0.14 w of −0.43 w), and the band adds ±(¼ of it + 0.02 w) |
| depth weight between the two limits | G | taken from the exact curved-duct solution of a reference bend r_m = w as the fraction of the way from its deep to its shallow conductance at the device's h/w (h/w 5: 0.05; 1: 0.33; 0.5: 0.58). It is a proxy rather than a 3-D result, and the record's band moves it by ±max(weight, 0.1). |
| square-corner term c_corner (rc7, square quarter turns only): −0.24 w (sharp) to −0.56 w (inner rounding r_in/w 0.36), deep and shallow limits blended by the depth weight; the outer rounding has no measurable effect to r_out/w 0.307 | G | a 2-D computation (`provenance/rc7_scripts/corner/`, Hele-Shaw and planar Stokes, with the solver checked on rc6's sharp-corner fits and exact arc), relative to the centerline through the corner point (so negative, unlike a U-turn fold's L_eq). It is the glossary's bend-equivalent length L_eq (EXTERNAL, supplied by this work, G:1448), computed rather than calibrated (spec v5 S:409), so no validator in the record discharges it. The record's band at the square-turn head gives 121 rows beyond ±1%, 4.3% of F (LIMITATIONS 8). |
| 1-D Kirchhoff network on the drawn widths | none | the model rather than a measurement, so 3-D effects at gap mouths and the outlet bifurcation are outside it; the junction terms above are fits whose range the record reports (model band), and their 3-D calibration is future work (docs/LIMITATIONS.md) |

## Gifford 2014, Fig. 4 bead device (`gifford2014_fig4_design`, `gifford2014_fig4_as_fabricated_{20,8}um`)

| value | code | source / note |
|---|---|---|
| D 5.2032 µm | D | printed f_gap 3.34 × 10⁻⁴ at w_c 100 µm (P) through κ_g |
| F 0.90644 | D | the outlet identity at the printed final side width 300 µm; the march overshoots it by one row (0.906474) |
| depth 80 µm, w_c 100 µm constant | P | Fig. 4 caption (80 ± 3 µm) |
| G 20.5, g_p 25.3 µm | M | G is the record's measurement on Fig. 4, and g_p = the inset's row pitch minus G, with the pitch 45.82 µm (four post columns, 45.81–45.84) at its 300 µm bar read with its anti-aliased ends (107.52 px), rc6 (`provenance/rc6_scripts/measure_rc6.py`). In rc5, the bar was read at its full-intensity pixels only (107 px), which gave a pitch of 46.1 and g_p 25.6. |
| l_g 59.5 µm (post across the flow) | M | l_g is the separation of each leg's two post columns on the Fig. 4 inset, 159.5 µm center to center (159.8 and 159.2 on legs 2 and 4), minus the printed 100 µm center (by dark-line centers across, 59.7 ± 0.3 µm on 96 posts; leg 4's columns alone, 59.2 µm; 59.5 is the value used). The center-to-center rule is free of the outline bias that relief contrast gives (rc6). In rc5, the reading was 56.8 between dark outline centers. |
| 4 legs, nested (spiral-out) paperclip, three ports at the top | P + M | caption "the second and fourth legs of the device"; schematic (not drawn to scale) |
| leg spacing 818 µm (all three Fig. 4 presets) | M | In rc6, the center lines of legs 2 and 4 on the Fig. 4 inset (neighbors in a 4-leg opencif paperclip) were found 817.9 µm apart at its bar (`provenance/rc6_traces/code/panels_gifford.py`), whereas the generator's rule gave 969 µm. On the rebuilt mask, leg 2's posts and center line sit on the image's, but its side channels there read about 18 µm wider per side than the rows that the mask places at that point (`measured.leg2_side_walls`; the nested layout fixes which rows pass the inset, 1,730 to 1,758; register GIF4-1, open). |
| leg 4's array ends 2 mm early (`pc_last_short`) | M | inset reading; *reproduction-only* |
| port radii 750 / 300 / 300 µm | A | the schematic is not to scale |
| sealed start, w_s(0) = 0 (design-level) | P | p. 4498 and Fig. 2 convention. The start draws a 4.86 µm first side channel (the pointed side-channel starts of the zooms' panels B and C), narrower than the particles that the first gaps are meant to pass (up to D = 5.20 µm), so the design-level drawing shows the equations' geometry rather than a device to make. The as-fabricated presets start at 20 or 8 µm. |
| floor 4.5 µm, corner limit 1.8 µm (design-level) | G | chosen after the design was known, so that the 4.86 µm channel would pass, and therefore not a fabrication fact |
| start side width 20 µm / 8 µm (as fabricated) | A | Gifford prints neither w_s(1) nor a floor. A side width of 20 µm is the student's floor (his gaps are about 20.5 µm), and 8 µm is the "10 : 1 aspect ratio limitation" that he mentions when discussing pore size against depth (p. 4503), applied to 80 µm, which is an inference rather than a stated floor. |
| start rule (as fabricated) | P | p. 4498, which says to begin with w_s(1) "close to or equal to zero, depending on the device production capabilities available", recurse eqn (4), and stop at the final width; modeled as a start side flow S0 in the design and a sealed drawing, stopping on the delivered F |
| flow rate (Reynolds gate) | A | not printed for Fig. 4; default 30 µL/min |
| no meander | P | Gifford has none; meander segments are Xia 2016's addition |

In RC8, `gifford2014_fig4_aligned` is the nested reading (register GIF4-1 / INV-2; D1 of the adjudication of the rc7
reviews). It keeps every printed value and takes the layout that the caption, the schematic, Fig. 6, and the inset's
taper direction support (and the inset's arrows, the weakest line), and the frozen physics is exact on the drawing,
with a declared off-design entrance. It does not reproduce leg 2's side channels beside the inset, which are 17.5–18.1
µm per side wider on the image than the printed schedule draws them (open). It differs from the design preset in the
following values.

| value | code | source / note |
|---|---|---|
| start 8 µm (`start: floor`) | A | Gifford's rule (p. 4498) at his stated "10 : 1 aspect ratio limitation" (p. 4503) at 80 µm, as in `gifford2014_fig4_as_fabricated_8um`. This is App. A's "other start mode" (A:1118–1119; spec v5 S:271), a fictitious start flow 2σ0 = 0.0022 of the inlet's, realized on the sealed drawing as rows 1–3 off-design (row 1 2.0 D, 0.5% of F; declared). The 8 µm row-1 side channel is narrower than the 10.2 µm particles that rows 1–3 pass (and than Gifford's 8.3 µm beads, of which the streamline model sends 0.1% there). The inset's beads (INV-2) allow at most about 5% early loss and bracket a straight start at up to about 12 µm, so 17.2 µm (the inset's own gap) and the lab's later 19–20 µm floor are rejected (LIMITATIONS 28). A 20 µm floor with the narrow rows folded into a head is the variant `provenance/rc7_traces/variants/gifford2014_fig4_fold20.yaml`, which is a figure rather than a reproduction. |
| floor 8 µm, fillet 2.5 µm | A / G | the floor at the start width (the 10 : 1 limitation at 80 µm), and the fillet below half of it |
| G 17.2, g_p 28.6 µm (row pitch 45.8 µm unchanged) | M | the inset's posts' flanks along the flow, read at search half-windows of 3.5 and 6 px (9.8 and 16.7 µm, both shorter than the gap) and checked at the trace window (3 px, 8.37 µm, at most half the gap; register GIF4-2, GIF4-3) |

The layout values are the design preset's, with four nested legs 818 µm apart (code M) and leg 4's array 2 mm short
(code M, the inset's leads). Traced on the inset (window 8.37 µm), the claim (leg 4) gives a median of 0.61 µm, a p90
of 1.40 µm, a coverage of 0.999, censored null ratios of 8.66 / 6.46, and chamfer null ratios of 2.89 / 2.38. Four
lines of evidence support the nested layout. The first is the caption ("the second and fourth legs"). The second is
the schematic's topology, because its "not drawn to scale" (p. 4502) is about scale rather than about which leg
neighbors which. The third is the inset's two arrows, which point the same way (one per visible channel if they follow
the same paper's Fig. 6, which draws one arrow per leg); their placement, symmetric about leg 4, also allows both to
mark leg 4, so this line is the weakest. The fourth is the taper direction (V13,
`provenance/rc7_traces/code/taper_sign.py`), by which the neighboring leg's outer width on the image narrows in the
direction of the nested leg 2's (PASS, with the image's slope −0.00240 ± 0.00059 µm/µm (4.1 SE) against the mask's
−0.00257, 130 of 130 eligible positions at a ±30 µm window) and against reading B's leg 3 (FAIL, −0.00254 ± 0.00062
(4.1 SE) against +0.00301, 117 of 117), with leg 4's own array part agreeing (PASS, −0.00469 ± 0.00040 (11.6 SE)
against −0.00800, 80 of 80 at ±12 µm). The reading does not reproduce leg 2's side channels. Drawn at 105.1–106.7 µm
on the rows that pass the inset (1,730 to 1,758), they read 122.6–124.8 µm on the image (the drawn widths plus the
walls' measured offsets; the width profiles read 124.5 / 125.0 µm), which is what the schedule draws about 280 rows
(12.7–13.0 mm) further on (1.17 ×). No reading with the printed single f and a start inside the beads' bracket
explains it.

Since RC8, rc7's `gifford2014_fig4_aligned`, reading B, is the variant
`provenance/rc7_traces/variants/gifford2014_fig4_pathorder.yaml`, a rejected reading that is built, traced, and shown.
It places the legs side by side in path order (`layout.order: path`, code A, and P contradicted by the caption, the
schematic, and Fig. 6's nesting), so that the inset's neighbor is leg 3, whose side channels (124.7 µm, code M) the
printed design reaches when leg 4's array is made 8.2 mm short (`pc_last_short` 8200, D, fitted to that width). Traced
on the inset, the claim gives a median of 0.62 µm and a p90 of 1.40 µm, and the neighbor's mean wall offsets lie within
1.4 µm, which is consistency rather than a test (because of the fitted end) and a mean that cannot see a slope. Its
taper runs against the image's (V13 FAIL). The rc6 Fig. 4 package's device (a 51 µm start, nested) fits the same
widths only through a start flow that the inset's beads rule out (INV-2), and it is kept as
`provenance/rc7_traces/variants/gifford2014_fig4_start51.yaml`. The nested reading with the design preset's sealed
start (every row at D) would draw a 4.86 µm first side channel, 16.5 : 1 at 80 µm depth, against the printed 10 : 1
limitation, whereas the 8 µm start keeps the printed statement at the cost of three declared rows.

## Gifford 2014, Fig. 5 (`gifford2014_fig5{A,B,C,D}_as_fabricated_12p5um`)

| value | code | source / note |
|---|---|---|
| D 4.2200 / 5.3268 / 6.9823 / 9.0268 µm | D | printed f_gap 2.2 / 3.5 / 6.0 / 10.0 × 10⁻⁴ at w_c 100 µm through κ_g |
| F 0.95748 | D | the outlet identity at the printed final side width 450 µm (text "~23 : 1", caption "24.4×") |
| depth 125 µm, w_c 100 µm constant, 25 µL/min | P | Fig. 5 caption and text |
| G 20.5, g_p 25.3, l_g 59.5 µm | A | not printed for Fig. 5; taken from Fig. 4 (same lab and period; rc6 values, rc5 25.6 / 56.8) |
| layout, a nested paperclip fitted to 30 × 50 mm | A | Fig. 5's devices are not drawn |
| start side width 12.5 µm | A | the 10 : 1 remark (p. 4503) applied to 125 µm |
| port radii | A | |
| measured comparison | P | Fig. 5 enrichment by bead size (Coulter counts, caption). This project's reading (`provenance/gifford_fig5_d50.py`; three interpolants × three full-retention ceilings, which are the duct law's 23.5×, the caption's 24.4×, and the text's ~23 : 1) gives the retained fraction at D, R(D) = 0.39–0.42 / 0.44–0.54 / 0.38–0.48 / 0.17–0.30, and d₅₀/D = 1.09–1.14 / 0.98–1.03 / 1.01–1.04 / ≥ 1.10 (5D above 10 µm under the caption's ceiling). The earlier, unshipped figures 1.135 / 1.006 / 1.027 / ~1.108 are this reading at the caption's (5A to 5C) and text's (5D) ceilings. The duct-law ceiling gives R > 1 for the 10 µm beads of 5A and 5B. |

## Gifford 2014, Fig. 6 platelet device (`gifford2014_fig6_design`, `gifford2014_fig6_as_fabricated_15um`)

| value | code | source / note |
|---|---|---|
| D 3.6236 µm | D | printed f_gap 1.04 × 10⁻⁴ at w_c 125 µm through κ_g (u = 0.029, inside the κ knee) |
| F 0.71768 | D | the outlet identity at the caption's final side width 140 µm (the text says 144, which gives F 0.72923) |
| depth 150 µm, w_c 125 µm, 55 µL/min, ~40 cm, < 6 cm² | P | Fig. 6 caption, pp. 4501–4502 |
| 7 nested legs, inlet inside the innermost loop, ring-like outlet | M | Fig. 6 schematic, which Gifford says is not drawn to scale |
| leg axes `reproduction.leg_Y` [0, 3491.4, −587.8, 4079.3, −1175.4, 4666.8, −1764.3] µm | M + A | In rc6 (both Fig. 6 presets), legs 1, 3, 5, and 7 lie on the micrograph's post rows (588 µm apart; M), and legs 2, 4, and 6 lie at 2 × 1,745.7 µm − Y of the leg before, where semicircular turns about the 16 turn walls' common center put them (the center is M, with a residual SD of 0.13 px and radii of 1.68–3.17 mm, `provenance/rc6_traces/code/gifford6_turns.py`; the semicircle is A, opencif's turn model, because only 8–23° of each turn is in the frame). The generator's uniform rule (723 µm, turns of radius 362 µm) was far tighter than the image's turns; *reproduction-only*. |
| last leg's array 257 µm short (`reproduction.pc_last_short`) | M | In rc6, the micrograph's leg 7 ends 5 rows before legs 5, 3, and 1, and the generator ends it 4 rows beyond the image's (it moves the end in steps of about n − 1 rows, LIMITATIONS 25); *reproduction-only*. |
| G 20.5 µm | A | as Fig. 4 (the micrograph reads 17.1 µm along the flow, where relief contrast shifts outlines, and the reading is not used) |
| g_p 30.9, l_g 70.8 µm | M + A | Fig. 6's micrograph at its printed 125 µm center, read between the inner tip lines of each leg's two post rows (three legs; synthetic calibration of the tip estimator), gives a row pitch of 51.38 µm (51.29–51.41; at most 52.01 by the caption's 250 µm bar, which runs off the raster) and l_g 70.8 µm (70.6–70.9), and g_p = the pitch minus the assumed G (rc6, `provenance/rc6_scripts/measure_gifford_fig6.py`). In rc5, l_g was 56.8 as for Fig. 4, and g_p was 32.2 from an unscripted pitch reading of 52.7 µm, above the bar's bound. |
| ring and port sizes | A | parametric ring |
| floor 4.5 µm, corner 1.8 µm (design-level) | G | chosen so that the 4.69 µm sealed first channel draws |
| start side width 15 µm (as fabricated) | A | the 10 : 1 remark applied to 150 µm |
| measured comparison | P | platelet recovery 83.2% (80–85%; Gifford's "yield", which is not App. B's F), ~3× enrichment, side : center ~2.5× |

The rc7 preset `gifford2014_fig6_aligned` (the rc6 Fig. 6 package's device, version 2) differs from the design preset
in the following values.

| value | code | source / note |
|---|---|---|
| F 0.72923 | D | the outlet identity at the text's 144 µm (p. 4502; the caption's 140 µm gives 0.71768), and the micrograph reads 145.5 µm at leg 7's end. The paper's platelet counts give F ≈ 0.69 in operation (register INV-3). |
| floor start, w_s(1) 40 µm | D (fitted) | the start that puts leg 1's side channels at the micrograph's 58.9 µm (code M). Legs 3 and 5 test it in part (the text's 144 µm, which sets leg 7's end through F, was preferred over the caption's 140 µm with legs 3 and 5 in view, GIF6-2; leg 7 is not a test). Its declared entrance is rows 1–40 above D, with row 1 at 9.9 D (35.9 µm) and 14.5% of F. Its early side flow (2σ0 = 0.11 of the inlet's) is consistent with the 83% platelet recovery (Gifford's "yield"; App. B reserves 'yield' for F, B:4199–4201; register INV-2). |
| G 17.1, g_p 34.3 µm (row pitch 51.4 µm unchanged) | M | the micrograph's posts along the flow (rc6 `gifford_fig6_measurement.json`), used here rather than only reported |

## Dinh 2024 V1 (`dinh2024_v1_specified`, `dinh2024_v1_depicted`)

| value | code | source / note |
|---|---|---|
| D 2.1742 µm | D | printed f* = 1.04 × 10⁻⁴ at the inlet width 75 µm through κ_g |
| F 0.5592 | D | the frozen outlet identity at the printed outlet widths; measured 0.505 ± 0.033 (Fig. 4) |
| depth 140, G 19, GS 28, w_c 75 → 52, w_S(in) 21 µm, width 3.2 mm | P | Dinh p. 916 and Fig. 2 (p. 917) |
| outlet side width w_S(out) 44 µm (V1; 48 and 57 µm for V2 and V3) | P | Dinh 2024; with w_c 52 µm it sets the F target 0.5592 through the outlet identity |
| taper form, constant maximum shear between 75 and 52 µm | A | Dinh prints only the end widths |
| g_p 26.1, l_g 54 µm | M | The Fig. 2B micrographs, at each panel's 50 µm bar read with its anti-aliased ends (46.91 px; rc5 read 46 px), give the row pitch G + g_p = 45.11 µm (six post-column readings on panels 2–4, 45.05–45.17, SD 0.04; ± 0.2 µm with the bar's ± 0.3%), from which the printed G 19 is subtracted. l_g is the post columns' separation, 106.06 µm center to center, minus the printed 52 µm (rc6, `provenance/rc6_scripts/measure_rc6.py`). In rc5, g_p was 27, from a pitch of 45.98 ± 1.0 µm (`provenance/rc5_scripts/measure_dinh_pitch_gap.py`). |
| meander wall 27 µm, three goal-post arms at R0 1,281 µm | M | Fig. 2A inset |
| square quarter turns, outer corner rounded 7 µm, inner 10 µm (depicted; rc7) | M | On Fig. 2B's first and second panels, the head's quarter turns are right-angle corners rather than rc6's 41.5 µm-radius arcs (register DINH-1). They are drawn with `meander.turns: square` and the corner term c_corner from 2-D solves (`provenance/rc7_scripts/corner/`). The specified preset keeps the arcs (the text prints no corner shape). |
| two-stage inlet filter, port radii 1,020 / 980 / 956 µm | M | the Fig. 2A main panel, read at the arm span taken as the 3.2 mm chip width, which is a scale that is not established (the panel read at it is 60.6 mm long against the printed 75.8 mm, and the round port shows only that the panel is not stretched one way), so the radii carry that scale's uncertainty (rc5, after a review of rc4) |
| floor 19 µm | A | Dinh prints none; set equal to the printed gap, following Xia 2016's rule that a meander-headed CIF device's minimum feature is its gap G (line 99) |
| gap-1 fraction 21 f* (specified) | G | the smallest that lets three arms fit 3.2 mm; *departure*, with a gap-1 cutoff of 10.0 µm against 2.17 |
| rectilinear transition, side channel 28 µm beside 3 wall blocks (`transition: paper`) | M | Fig. 2B second panel; *departure*, up to 13 × D at those rows |
| head fractions as depicted (depicted) | M | the Fig. 2A inset unit lengths (drawn at those lengths, with a resistance that includes the corner terms of the square quarter turns and the bend correction of the U-turn arcs; in rc7, after code review R2, re-derived for the square-turn head, which moves the fractions by +0.27 … +0.95%, row 4 −1.31%, against rc6's arc-turn derivation, and makes the planned head 219 rows; `provenance/rc7_scripts/regen_depicted_dinh_square.py`); *departure* |
| flow 30 µL/min | P | Dinh 2024 (operating flow, used by the Reynolds gate) |
| lateral margin 15 µm | G | |
| resulting length ~150 mm against the stated 75.8 mm | none | In rc7 (register INV-1, independently checked), Fig. 2A, drawn to scale, and Fig. 3A's photograph put the linear section at 23–29 mm (514–642 rows; at most about 36 mm), where these presets draw 2,482 rows (112 mm). Through the closure at the printed seam and outlet widths and App. A R13(a)'s invariant lamina, those row counts give a linear-section cutoff of 4.28–4.79 µm at F 0.5592 (4.04–4.52 µm at the measured 0.505). This is a two-cutoff reading of the figures, offered as a second preset (`dinh2024_v1_figure_aligned`, below) on the author's decision, whereas rc7 final fits its cutoff to the printed length instead (3.67 µm, below). The rc6 analysis finds a conditional conflict (in the rc5 wording, corrected after the two reviews of rc5 and independent checks of that verification, with rc6 numbers at the re-measured pitch from rc5's scripts ported to take the preset's own pitch, `provenance/rc6_scripts/dinh_length_bound.py`, `dinh_taper_alternatives.py`, and `dinh_schedule_at_F.py`, which reproduce rc5's outputs exactly at 46.0 µm). There are four premises. The first is the frozen cutoff law's cubic kernel with a lamina held the same at every gap at the plateau prefactor κ_g, which is App. B's own reading of Dinh's rule ("Under the Xia-scaled reading, which is the one Dinh describes, the lamina is the same at every gap", App. B :2970–2972; Dinh adjusts each gap "to maintain the desired size cut-off throughout the device", p. 916), fixed by the printed f* = 1.04 × 10⁻⁴ at 75 µm (D 2.1742 µm on the κ plateau). At Dinh's rows wider than 62 µm, u = D/w_c < 0.035, so the law is below its plateau and only the cutoff's value is uncertain. Under this reading, App. B's band "collapses to" [1.80, 2.17] µm (App. B :2972; the records' 1.805–2.174 µm), whereas App. B's recorded band, which holds f at its initial value instead, is [1.25, 2.17] µm. κ drops out of the row count, and the 1-D duct law enters only the presets' own schedules in (2). The second premise is every center width at or above the printed 52 µm, so that a row extracts at most f_max = 2.1597 × 10⁻⁴. The third is two gaps per row (the march's recursion P(i+1) = P(i)(1 − 2 f(i))). The fourth is every gap row at least the measured 45.1 ± 0.2 µm pitch long (rc5 used 46.0 ± 1.0 µm, with the bar read short). (1) The lower bound for any taper (every row at f_max) is N ≥ ln(1 − F) / ln(1 − 2 f_max). At the F that the presets close (0.5593, the outlet identity at the printed 44 / 52 µm, a model value that is not printed), the bound is ≥ 1,897 rows, 85.6 mm (85.2–85.9). At F 0.5105 (the assumed constant-shear taper reaching 52 µm, code A), it is 74.6 mm. At the measured F 0.505 (Fig. 4 caption, p. 919, "50.5 ± 3.3%" and "Values shown are mean ± standard deviation (panel A: n = 3; panel B: n = 4)"; p. 916, "50.1 ± 1.6%" of the sample volume), it is 73.4 mm (73.1–73.8 over the pitch's ± 0.2 µm). Crediting the presets' reconstructed heads, whose first rows are drawn above f_max (the depicted head's first 18 rows, up to 6.4 f_max, and the specified head's gap 1, 10.1 f_max, a departure from the constant lamina), with N ≥ [ln(1 − F) − Σ over those rows of (ln(1 − 2 f_i) − ln(1 − 2 f_max))] / ln(1 − 2 f_max), gives 84.7 and 85.2 mm, 73.8 and 74.2 mm, and 72.6 and 73.0 mm, and ≥ 1,764 linear rows (79.6 mm) after the reconstructed head at 0.5593. These are the gap rows alone. The printed "∼3.2 mm × 75.8 mm" is "An assembled CIF device" (Fig. 3A caption, p. 918; also p. 916), which "included a single inlet with an integrated filter … and two distinct outlets" (p. 916), biopsy-punched (p. 914), of sizes that are not printed. At the measured F the bound leaves 2.4 mm for them (up to 3.2 mm with the heads' credit), and 1.2 mm (up to 2.0) at F 0.5105, so by itself it shows no conflict at either. (2) For a progressive taper, the bound is reached only with the center channel at 52 µm in every row. Dinh's linear section narrows it "progressively" (p. 916), "from 75 μm near the inlet to 52 μm near the outlet" (Fig. 2 caption, p. 917), and such a taper needs more rows. At F 0.505, it needs 2,346 rows (105.8 mm) for a width linear in the row index (`dinh_taper_alternatives.py`), and the presets' own schedules (their heads and the assumed constant-shear taper, code A) need 125.5 mm (depicted; rc6 125.6, before the square turns) and 128.6 mm (specified) of head and rows, 99.9 mm of them linear rows (`dinh_schedule_at_F.py`, planned at the target 0.5592, 0.5105, 0.505, and 0.472; in rc7, `provenance/rc7_scripts/dinh_lengths/`; 137.6–140.6 mm at the presets' F, and 118.4–121.5 mm at 0.472, one SD below the measured mean). The rows fit the printed length at F 0.505 only if their constant-equivalent extraction is at least 96.9% of f_max with nothing else in that length, which a constant 52.8 µm center would give (the condition is on extraction rather than on mean width, because f is convex in w). A linear narrowing followed by 52 µm fits only if it ends within the first 10.2% of the rows (11.7% at a 44.9 µm pitch, 8.8% at 45.3 µm; 11.9–13.6% crediting the heads; at F 0.472, 39.2%). Thus, under these premises the conflict holds at the presets' F for any taper, and at the measured mean F for a progressive taper at the measured pitch. It is smaller at the measured mean than at the presets' F (and smaller again at one SD below the mean). It is a conflict between the premised model and the printed footprint rather than among the printed numbers, and its cause is not known. The cause may be the kernel's lamina-flux relation in Dinh's 140 µm deep channels, the constant lamina, the pitch measurement, or a taper or design rule that the paper does not print. The outlet identity's F does not explain it, because at the measured F the bound's conflict goes whereas the taper's remains. Nor does the cutoff law below its plateau, because readings that hold D itself constant with the frozen Hill curve κ(u) below the plateau either make the device longer or do not return the printed f* (`dinh_taper_alternatives.out`, "below-plateau" lines). A D that returns f* at 75 µm with that curve is smaller (1.02 or 1.16 µm), and as the channel narrows κ(u) rises, so the lamina thins downstream, and the every-row bound alone is 132.1–142.3 mm at the measured F. The two readings that do not return f* (they give 1.27 × 10⁻⁴ and 1.51 × 10⁻⁴ at 75 µm) move the lengths both ways. D 1.805 µm, the lower end of App. B's band for a constant lamina, lengthens the bound to 85.6 mm and shortens the linear taper to 105.35 mm, whereas D held at 2.1742 µm leaves the bound unchanged and shortens the linear taper to 91.6 mm, still above 75.8 mm. Neither removes the conflict at the measured pitch, because a linear narrowing fits within none, or the first 15.7%, of the rows (18.0% at a 44.9 µm pitch for the second). In rc4, the documents stated "at least 1762 linear rows for any taper". Its p = 0.4672 was rounded (the exact recursion with the exact p gives 1,764), and it left F unstated as a premise (rc6, `provenance/rc6_scripts/dinh_length_bound.py` and `dinh_length_bound_{depicted,specified}.json`; rc5's at 46.0 µm in `provenance/rc5_scripts/`). |
| outlet end as fabricated. Each side channel's outer wall runs straight for 65.4 ± 2 µm from the last gap's center (64.9 and 65.9 µm; 74.9–75.6 µm from the last pill; rc6, at the 50 µm bar read with its anti-aliased ends; rc5 read 66.2 and 67.2 µm at 46 px) and then flares at 14.50° (left) and 14.43° (right) per side, each from its own wall's straight direction (24 fit intervals per wall, 14.46–14.56° and 14.30–14.56°), whereas the inner wall is straight. The wall direction from which the angle is taken is the least certain input. Fitted over six row ranges, the straight walls tilt 0.28–0.40° (left) and 0.51–0.73° (right), giving 14.43–14.55° and 14.38–14.60° per side. A review of rc5 measured the right wall at 0.09–0.35°, 14.7–15.0°, on an image at four times the native resolution, which is not reproduced on the native raster, where the "52 µm" label's glyphs stay at least 5 px from the wall's outline and the rows where they come closest are excluded. The center channel after the last gap is 49.7 ± 0.2 µm outline to outline (52.8 µm between the post columns, where the caption prints 52 µm) and begins to widen 191 µm after the last gap's center, at 0.25 µm per µm, to 58.5–58.9 µm at the panel's edge 228 µm after it (rc6 scale; rc5 read 50.7, 53.8, 195, 59.6–60.1, and 232 µm at 46 px). It is not drawn, because the center lead is the generator's (A). | M | Fig. 2B fourth panel, native raster (sha256 a21e2e3ce0710e32410ba4941f348c6a1ef561b748f569e24decce1dea9fbbe6; 1.066 µm/px from the 50 µm bar read with its anti-aliased ends, 46.91 px, about ± 0.3%; rc5 1.087 at 46 px; the two straight walls are tilted 0.35° and 0.70° in the raster over rows 1,590–1,750, so angles from the image's vertical, as rc4's 14.0° and 15.2°, carry that tilt with opposite signs); `provenance/rc5_scripts/measure_dinh_flare.py` (v2, with the right wall's rows beside the "52 µm" label excluded, because two picks of v1 had landed on it); for the tilt ranges and the center channel, `measure_dinh_centre.py` (after the rc5 reviews; every pick in `dinh_centre_measurement.json`, ±0.5 px per edge) |
| `outlet.flare_deg: 14.5`, `outlet.flare_start_um: 65.4` | M | the row above (rc5, 67); drawn since rc5 (25 µm tangent arcs at both kinks) and checked on the drawn mask by `tests/test_figure_anchors.py` |
| the lead width the flare reaches, 105.75 µm, and the teardrop lead after it | A | not visible in the panel; the generator's closure (the teardrop's compact-run rule) chooses it |

The rc7 preset `dinh2024_v1_figure_aligned` (register INV-1, a two-cutoff reading of Fig. 2A's head and the printed
75.8 mm length, whereas the printed `dinh2024_v1_depicted` stays the reference) differs from the depicted preset in the
following values.

| value | code | source / note |
|---|---|---|
| D 3.67 µm (the device's cutoff, set by its linear section) | D (fitted) | In rc7 final (the decision of A. H. Routt on 2026-10-19 to fit the printed length; register INV-1, revision), D is the cutoff at which the drawn device (the depicted head, the linear section, the inlet with its filter, and the outlets) is the printed 75.8 mm long, with 75.71 mm drawn and 874 linear rows (39.5 mm), the fit moving about 0.09 µm per 2 mm of length (`provenance/rc7_scripts/postfreeze/length_vs_D.py`). The closure at the printed seam (21 µm) and outlet (44 / 52 µm) widths fixes the linear section's total extraction with the mass identity (A:1126–1132, A:1141, A:1148), and given that extraction, App. A R13(a)'s invariant lamina (A:1294–1303; the rule that Dinh states, B:2958–2962) with the assumed taper (code A) turns the row count into the cutoff. Dinh's bead data (ESI Fig. S1B; App. B P-C7.2, B:2986–2996) put the cutoff below 4.7 µm. Against Gifford's measured soft edge (Fig. 5, four devices of known D), the 4.7 µm beads' passage (3.6, 6.7, and 14.2% of the fluid's at the three concentrations) implies D 2.3–3.9 µm, inside which 3.67 µm lies (Gifford's curves predict 11–35% passage at 3.67 µm, 29–50% at 4.29 µm, and 0–4% at 2.17 µm; `provenance/rc7_scripts/postfreeze/softedge.py`). The figures' row counts give a larger cutoff. Fig. 2A (drawn to scale) and Fig. 3A's photograph put the linear section at 514 rows (23.2 mm; Fig. 2A read at the printed 3.2 mm width, as in rc6's `data/r4/dinh_layout_meas.json`) to 642 rows (28.9 mm; at the printed 75.8 mm length), which gives 4.79–4.28 µm at the 1-D outlet identity's F 0.5592 and 4.52–4.04 µm at the measured F 0.505 (A:816–821). The first reading shipped in rc7 took 4.29 µm (642 rows, with 65.2 mm drawn, 10.6 mm short of the printed length). The refit gives up the figures' proportion for the printed length, because the drawn linear section starts at 41.5% of the inlet port → middle outlet distance, against 54.2% on Fig. 2A and 46–58% on Fig. 3A's photograph (whose bound puts the linear section at most about 36 mm, about 800 rows, where this preset draws 874), and Fig. 2A's scale-free head : linear proportion is 1.09 where this preset draws 0.65. No single reading of the figures gives the printed length, the depicted head, and the figures' proportion together (register INV-1). The linear rows reach f_gap 6.13 × 10⁻⁴ at 52 µm, above κ_g's calibration cap 5.76 × 10⁻⁴ on the 133 linear rows narrower than about 53.6 µm (21% of the linear section's extraction) and above spec v5's f_soft on 252 (`validation.gates.f_beyond_calibration`). A prefactor 10% off beyond the cap moves D by −0.05 to +0.02 µm (3.63–3.69 µm; `provenance/rc7_scripts/postfreeze/kappa_beyond_cap.py`). D is a cutoff, row by row, rather than an over-extraction factor (A:809–811). |
| the head: 219 rows prescribed (`reproduction.f_head`, `NM_force` 219; `provenance/rc7_scripts/make_dinh_figure_aligned_head.py`) | M + A | exactly the depicted preset's planned head, which is its 18 rows measured on Fig. 2A's inset (2.8–13 f*, mostly about 3.7 f*, cutoffs about 4.2 µm), then 201 rows at the printed f* (2.17 µm), which give the head nearly Fig. 2A's length at its lateral reading (24.1 mm drawn from gap row 1, against 25.4 mm on the figure, which is 1.3 mm or 5% longer; at Fig. 2A's length reading the figure's head is 31.8 mm, as in the row above); declared departures from the device's D ('head drawn as depicted'), below the κ plateau |

The preset contradicts Dinh's "to maintain the desired size cut-off throughout the device" (p. 916) and the printed f*
in the linear section (App. B reads Dinh's design as one lamina). Beyond the layout, it has support, which falls short
of proof, in the measured head units and in the filtrate's platelets (mean 6.8 fL; 4.9% of the 12 fL ones), which a
2.17 µm cutoff would not pass. Platelets are outside the frozen rigid-sphere scope (B:4194).

## Xia 2016 RS-7 (`xia2016_rs7_specified`, `xia2016_rs7_depicted`, `xia2016_rs7_depicted_pullR`)

| value | code | source / note |
|---|---|---|
| D 6.8644 µm | D | printed f_ref 5.8 × 10⁻⁴ at w_c(ref) 100 µm (corrected eqn 1, 2022 correction) through κ_g |
| F 0.938 | D | the printed ~15 : 1 filtrate : retentate. Fig. 3g's last side channel reads 624.6 µm at 1.798 µm/px, against the 600 µm that this F draws, and F 0.94072 (15.9 : 1) with D 6.5705 µm (Fig. 3a's legs filled) gives both. That pair is the rc6 variant `provenance/rc6_traces/variants/xia2016_rs7_depicted_fig3.yaml` (`code/xia_fig3_fit.py`), a reading of the figures rather than of the text. |
| G 16, GS 22.8, w_c 300 → 150 µm at constant maximum shear, w_s(1) = G | P | Xia pp. 3–4, Fig. 1. The principle and the ends are printed, whereas the profile between them is not, so opencif draws its constant-shear form (A in its form; no user profile, LIMITATIONS 16; the panels' widths, LIMITATIONS 24 c). |
| depth 140 µm | A | printed 135–150 µm |
| g_p 35.8, l_g 70.9 µm | M | Figs. 3g and 3c at 1.798 µm/px, which is the scale that their printed widths (w_c 300 at the first posts, 150 at the outlet), their post columns, and the curvature of panel d's loop-4 axis arc (1.7943 µm/px; the seven-arc fit of `provenance/rc6_traces/measurements/xia_fig3d_arc_scale.json`, 1.7922 ± 0.0044, would give a mean of 1.797) give (the single 250 µm bar, drawn on panel i, gives 1.936). g_p is the row pitch, 28.80 px = 51.78 µm, minus G 16, and l_g is the post columns' center-to-center separation minus the printed center (70.89 at the outlet, 71.00 at the first posts) (rc6, `provenance/rc6_scripts/measure_rc6.py`). In rc5, the values were 40.25 / 78.1 µm, at panel i's bar read as 1.953 µm/px. |
| meander wall 35.8 µm, head width 1.84 mm, hairpin (serpentine) head | M | Fig. 3b's head unit pitch and Fig. 3c's transition block pitch, both 65.15 px = 117.1 µm = 2 (GS + wall) with the printed GS 22.8; the head envelope at 1.798 µm/px (rc5, 41 µm and 2.0 mm); the serpentine start is Xia's own |
| 7 leg axes, inlet and outlet ports, ring outlet geometry | M | Fig. 3a at its 5 mm bar read with its anti-aliased ends (178.32 px, whereas rc5 used 179, so that every Fig. 3a length is scaled by 1.003827) |
| nine rectilinear rows, side 22.8 µm beside the blocks | M | Fig. 1c(i, ii); *departure* (Xia's deliberate pull into the center) |
| loop pull-in 2 µm | A | Xia writes "subtract 1–5 µm" (p. 4); Fig. 1e(ii)'s printed loop-1 widths point to about 1.5–2 µm (57.6 + 77.4 = 135.0 µm against 2 × 69.5; 2.0 under a linear rule, 1.5 with the model's own loop-1 offset of 1.0 µm), a consistency inference; *departure* |
| outlet tubing I.D.s 0.25 / 0.86 mm | P | p. 10 |
| flow 800 µL/min | P (range) + A (value) | Xia ran RS-7 at "~1 mL/min" (line 18) and ">0.75 mL/min" (line 244), and 800 is our choice inside the printed range. At this flow Re exceeds the gate of 10, so the device is outside the validity envelope on Reynolds as well as on the κ plateau. |
| tubing length 0.15 m | A | not printed |
| floor 16 µm | P | Xia 2016 line 99, "the minimum feature size ... is now equivalent to its gap size, G" (and line 62, "~20 µm" in general) |
| corner limit 4.5 µm | G | |
| gap-1 fraction 3.16 f* (specified) | G | the smallest that fits the 1.84 mm head; *departure* |
| head as depicted (depicted) | M + A | Figs. 3b and 3a at 1.798 µm/px; 10 units measured, 19 extrapolated, 76 rectilinear rows assumed (rc5, 21 extrapolated and 65 rows); drawn at the depicted lengths, with the resistance including the bend correction; *departure* |
| `loop_pull_mode: resistance` (pullR) | G | an option-combination check rather than a published device |

The rc7 preset `xia2016_rs7_aligned` (rc6's traced variant `xia2016_rs7_depicted_fig3` with the parallel-plate taper
that Fig. 3's panels d and e read, a reproduce-mode departure from spec v5 S16-b) differs from the depicted preset in
the following values.

| value | code | source / note |
|---|---|---|
| `centre_channel.shear_formula: parallel_plates` | A (chosen with panels d and e in view, which then give consistency rather than a test) | Xia holds "a relatively constant maximum shear" along each leg (pp. 3–4), computed with "formulas for rectangular channels based on Yang et al."; the glossary identifies Yang's Eqn 3 as the rectangular-duct resistance series (G:2036–2038), and Xia's text does not state the shear coefficient. Spec v5 holds constant shear with C(α) (S:456) and reports it "via C(alpha) [NEVER 6]" (S:606). Panels d and e's leg-3 center widths (219.8, 186.0 µm) lie on the parallel-plate form w_c ∝ Q_c (218.1, 186.1) and 11–20 µm off the exact one, which makes the form a reading of the figures and a reproduce-mode departure from spec v5 (register XIA-1). The record reports the true peak shear via C(α) (`validation.gates.shear_peak`), which rises 17.4% along the taper (`over_the_taper.max_over_min` 1.174) and is 1.34–1.57 times the parallel-plate value. |
| D 6.4159 µm, F 0.94072 | D (fitted) | refitted with that taper to Fig. 3a's legs and panel g's last side channel 624.6 µm (`provenance/rc7_scripts/xia_fig3_fit_pp.py`), as rc6's variant was fitted with the exact shear (D 6.5705 µm); the printed f_ref gives 6.8644 µm and ~15 : 1 gives F 0.938 |

## Meander head styles and where they come from

| style | origin |
|---|---|
| `hairpins` (serpentine segments, no arms) | Xia 2016 |
| `arms` (goal-post arms around the inlet) | Dinh 2024 Fig. 2A |
| `pack` (as many arms as the width allows, then hairpins) | this project |
| ring outlet | drawn after Xia 2016 Fig. 3a; Gifford's Fig. 6 schematic (not to scale) shows a similar ring-like outlet |
| `narrow` (units at the process floor where the full channel width would not fit) | this project |

## Examples renamed in this round

| file | why |
|---|---|
| `examples/wbc_paperclip_teardrop.yaml` (was `wbc_paperclip_ring`) | the ring outlet cannot carry this device's about 1.8 mm filtrate lead around 1 mm ports (it was refused), whereas the teardrop can |
