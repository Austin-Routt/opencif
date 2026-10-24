# Limitations and future work

We state everything below as a limitation of the frozen code base, and none of it is fixed there. The future work that
an item names needs 3-D computation, measurement, or new code that is outside the dissertation. Each item names what the
record already reports about it, so that a reader can see its size on any device.

## A. Physics scope (App. B thesis, Leg 5)

1. The physics is scoped to rigid, approximately spherical particles in a Newtonian, nearly constant-viscosity fluid.
   App. B's thesis says "rigid spheres in Stokes flow" (its Leg 5 reads
   "Rigid spheres in Stokes flow at nearly constant viscosity"), and the glossary states the regime by scale in the
   words "The framework is therefore not globally Stokes: it is scoped to rigid, approximately spherical particles in a Newtonian, nearly constant-viscosity fluid, with creeping-flow conditions at the gap and lamina scales that govern particle selection."
   The channel flow is laminar, although it is not creeping (the M14 channel gate allows Re ≤ 10, and the strict designs
   here run at 3.6 to 5.6). App. A explains why the 1-D model holds anyway. In a fully developed unidirectional channel,
   the Navier-Stokes equation "reduces exactly" to the Stokes equation
   "whatever the Reynolds number, as long as the flow is laminar and fully developed". Deformability, cell crowding, and
   whole blood are out of scope. Because platelets are discoid, platelet-specific validity is not established. D is a
   rigid-sphere diameter, and the equivalent diameter of a non-sphere is not established (App. B states that it is
   orientation-dependent). Examples named for platelets or white cells are rigid-sphere designs at those sizes.
2. The gates use one Newtonian fluid, which is water by default (ρ 1,000 kg/m³, μ 1 mPa·s). Since rc7 (adjudication item
   ADJ-C5), the `fluid` section (`density_kg_m3`, `viscosity_Pa_s`, `rheology: newtonian`) sets the Reynolds and Stokes
   gates. The gates scale with ρ/μ, so a more viscous fluid (plasma) runs lower. The schedule does not depend on the
   fluid, because its resistance ratios are independent of μ for a Newtonian fluid. No non-Newtonian model (blood's
   shear thinning, cell-free layers) is offered, and a value that the user gives is not checked against the fluid that
   it names. App. B notes that the channel gate binds for water at 90 µL/min and that, for platelet-rich plasma at that
   rate, the verdict turns on the viscosity. Give the viscosity in that case. The record reports the fluid in
   `validation.gates.fluid`.
3. The cutoff law's calibration. κ_g = 4.9138 [4.750, 5.078] (±3.4%) comes from Gifford 2014's array, with f_gap up to
   5.76 × 10⁻⁴ and u ≥ 0.035. Above that f, the law is extrapolated (the record gives `f_design_max`, with a note when
   it exceeds the cap). Below the plateau, the cutoff is a band by App. B C7's method, which is plateau-anchored rather
   than the implicit fixed point (App. B C7 CQ3) and which App. B warrants only as a low-severity consistency band
   (CQ4). The ±1% design tolerance is tighter than the law's own interval.
4. D computed on the exact cubic root lies at most 0.47% from κ_g's leading-order fit up to the calibration cap
   (tested).
5. F is the 1-D filtrate fraction. The 3-D flow at the outlet bifurcation (the glossary's η = 0.903, a measured
   correction) is not applied, so Dinh's network F is 10.7% above his measurement by construction (App. A R8; App. B
   P-C4.2). Since 1.0.0, every record reports ηF beside F, as the design code does (the record's
   `validation.design_code_reports.B6_F_real`), and the correction is reported and never applied.
6. CIF_v2's E11 and M16 are gates (rc3), and the physics auditor runs on every photomask (1.0.0), whereas its stricter
   P8 is reported rather than gated. Release rc3 ports CIF_v2's thin-lamina limit δ/h ≤ 0.10 (E11) and its
   quasi-developed-taper limit |w_c(i) − w_c(i−1)|/(G + g_p) ≤ 10⁻² (M16), with CIF_v2's own tolerances, into the
   validity envelope, where strict mode refuses and the other modes report (the record gives `lamina_over_depth` and
   `taper_rate_max`). Neither App. A nor App. B cites either limit. On the shipped designs, δ\*/h ≤ 0.016, and the taper
   rate is ≤ 1.34 × 10⁻³ inside the envelope and 6.4 × 10⁻³ (64% of the limit) on Xia RS-7 depicted, which lies outside
   the envelope on M13. CIF_v2 applies E11 only where it computes f from δ\* (its S1 path, so a prescribed schedule is
   skipped), whereas spec v5 writes it "∀i". Since rc4, the envelope takes spec v5's reading, in which E11 holds on
   every planned row, including the rows planned away from f\* (declared departures). The record gives
   `lamina_over_depth` on the cutoff rows, for parity with CIF_v2, and `lamina_over_depth_planned_max` with
   `lamina_all_rows_ok`. The planned maximum is ≤ 0.034 on the shipped designs, so no outcome moved. M16 has not been
   made to fail alone on a constructible design. On the shipped designs, the taper rate is 0.006 to 0.057 times M13, so
   M13 < 0.05 keeps it below about 3 × 10⁻³, and the seeded steep taper fails M13 first. The negative suite therefore
   seeds M16 alone on the gate function and, as an integration case, together with M13.

   Since 1.0.0, the physics auditor is ported and runs on every photomask. App. B gives the design code the channel gate
   Re_max = 10 (M14) and says that
   "the physics auditor separately checks the drawn device's inlet against its own, stricter channel limit" (P8). CIF_v2
   has both limits. Its `design()` gates Re ≤ 10, and its physics path requires an inlet Re ≤ 1.0 (a settable default,
   for which App. B gives no value) at the device's flow (`Q_in_uLmin`, default 5 µL/min). That requirement sits in the
   auditor's check `P8_stokes_gates`, which adds a naive entrance length 0.05 Re D_h ≤ Λ/2, and in the unified
   compiler's pre-flight Reynolds gate. CIF_v2's validators Z1 to Z3 and T1 require its demo devices, at 2 and 5 µL/min,
   to pass. Since 1.0.0, opencif runs that auditor, ported (`opencif/checks/auditor.py`, `docs/PHYSICS_AUDIT.md`), on
   the photomask file of every build, read back from the DXF. P0 (the file is the drawing), R0 (its own ruler on the
   straight rows, described below), P1 to P6b, P9, R6b (the posts), and P10 are pass conditions in every mode, and P7 is
   reported, as in CIF_v2. The design code's own gate M14 (Re ≤ 10 on the channel, and Re_gap and Re_δ ≤ 10⁻²) is
   evaluated on the drawn mask as the check M14, and outside it the verdict is exit 3. P8 itself is checked on the drawn
   inlet at the design flow and at 5 µL/min and recorded under its own name (the record's `validation.physics_audit`,
   check P8, with the planned values in `gates.P8`), although it is not a pass condition. At their own flows (25 to
   800 µL/min), every opencif design exceeds P8's Re limit (an inlet Re of 3.6 to 5.6 on the 27 built strict designs,
   6.7 on Gifford Fig. 6, and 61 on Xia RS-7), so a P8 gate would refuse them all. Among them are the Dinh V1 presets at
   30 µL/min (inlet Re 4.65), an operating point that App. B's scope discussion puts inside the design code's channel
   gate (a worst-case channel Re ≈ 5.2, against Re_max = 10). CIF_v2 runs its own lineage demo device at 2 µL/min. At
   5 µL/min, every design passes both criteria (an inlet Re of 0.38 to 0.93, and a naive L_e of at most 0.18 of Λ/2).
   Furthermore, the L_e criterion uses the formula that CIF_v2's own `entrance_length` calls "simply the wrong formula"
   in this regime, which spec v5 also rejects as a gate and which the glossary calls "wrong here". The port keeps every
   check's purpose and threshold. Where a DXF file shows something other than CIF_v2's pattern objects, the port uses
   the formula that the file supports. Thus, P1 anchors the cubic pair that opencif uses (whereas CIF_v2's physics
   module anchors its laminar kernel), P6b uses App. A's march-outlet identity (because CIF_v2's corrected closed form
   exists only for its canonical taper), and P9, R6b, and P10 are measured on the polygon (because a DXF has no strokes
   or fillet primitives). `tests/test_auditor.py` runs the port and CIF_v2's own auditor on CIF_v2's three demonstration
   devices, requires the same numbers for the network, the terminal window, P3 (its PE part), P4, P5, P7, and P8 and the
   same verdicts, and seeds fourteen faults, which exercise every hard and envelope check except P2 and P6b. App. B's
   sentence that the C7/C8 software is
   "warranted by their acceptance suite and the physics auditor" (lines 4375 to 4377) refers back to
   "the physics auditor" of lines 4349 and 4350, which App. B's provenance tag at lines 4371 and 4372 lists beside
   `CIF_v2.py`'s `reynolds` and `design` (M14). Up to RC8, opencif did not run it (the record's
   `assessment.status.physics_auditor_P8` read "NOT RUN ON THE MASK"), and open item D95 proposed naming CIF_v2 there.
   In 1.0.0, the sentence holds of opencif's auditor, with one difference that the dissertation should state, which is
   that the auditor's stricter limit is checked and reported whereas the design code's M14 is enforced. D95 is therefore
   closed into D106, revised for 1.0.0, which names opencif 1.0.0 and says which gate its auditor enforces
   (`docs/DISSERTATION_NOTES.md`, section 4; the frozen text is not edited here).

   The auditor has its own ruler. The generator measures its drawing where its plan puts each part and takes the parts'
   lengths from the plan (`opencif.checks.measure`), so an auditor that used only that measurement would share its
   errors. The auditor therefore measures every ordinary straight row again (`opencif.checks.ruler`, which shares no
   code with the generator's measurement), reading the posts found among the file's islands, the gap openings, the wall
   thickness, the row pitch, the center width at the post tips and at the corners that the tips fix at the gaps, and
   the side width to the outer wall, all from the file.
   The two measurements must agree within 10⁻³ (check R0), and the auditor's network uses the ruler's resistances on
   those rows. The ruler does not measure the head units, the U-turns, or the outlet, which stay on the generator's
   measurement, model terms, and plan-taken lengths (item 12).

   Since 1.0.0, CIF_v2's mandatory reports and acceptance validators are reported on every build. CIF_v2's `design()`
   attaches reports to every device (the total pressure drop, F_real = ηF, the transfer bias, the sensitivity map, and
   the wall porosity), and its acceptance run (`run_acceptance`) computes an operating-regime contract
   (`V_operating_regime`: Re ≤ 1.0 at 5 µL/min, 0.05 Re D_h ≤ Λ/2, wall shear ≤ 1,500 s⁻¹, and ΔP ≤ 200 kPa), a steric
   screen (`V_steric`: D/G < 0.5, D/h < 0.3, and δ\*/h < 0.2), and a Monte Carlo width tolerance (`V_fab_tolerance`: σ =
   0.25 µm on every center and side width, 24 draws, a yield equal to the share of draws with |F − F_target| ≤ 3 × 10⁻³,
   and ok when the yield is ≥ 0.95). None of them is a pass condition in CIF_v2 (`ok` = the gold devices, G14, and the
   J-suite). Since 1.0.0, opencif computes each of them by the design code's own formula on every build and records it,
   again without making it a pass condition (the record's `validation.design_code_reports`;
   `opencif/checks/design_code.py`; and `tests/test_design_code.py`, which checks each formula against CIF_v2's own
   function). On the shipped designs: {{DCR_LIM6}} `docs/CIF_V2_REQUIREMENTS.md` lists every gate, validator, and
   acceptance criterion of `CIF_v2.py` and its status in opencif.
7. The aspect ratio has been a gate in strict mode since rc4. CIF_v2's floor is w_min = max(F_min, h/AR_max)
   (`Process.w_min`, "the smallest buildable channel width"; spec v5, section 18). The glossary lists AR_max as an
   external fabrication limit, and Gifford 2014 assumes "a 10 : 1 aspect ratio limitation" (p. 4503).
   `process.aspect_ratio_max` defaults to CIF_v2's 10. Release rc3 ported w_min where CIF_v2 enforces it, in `design()`
   and its validator G7, which cover the gap G, the meander channel GS when a head is drawn, every center width w_c, and
   the terminal side width w_s(N). Spec v5 goes further, because G7 asks "widths ≥ w_min" of the drawing and G11 draws a
   side channel as a fold "iff w_s(i) < w_min". Because opencif draws no fold, rc4 makes G7 a gate on every straight
   side channel, both before drawing (`wmin_check`, on the planned w_s(i) of every straight row) and on the drawn mask
   (the check "aspect ratio of the narrowest CIF channel"). The drawn-mask check of rc4 took the gap openings at their
   design value and measured the straight and U-turn side channels. Since rc5, it takes the widths from the drawing,
   except the head units'. These widths are the narrowest measured gap opening of the linear rows, the center channel
   at both ends of its taper, every straight side section, the U-turn sides, and the three outlet leads
   (`validation.outlet_measured.lead_width_min_um`; a review of rc4 found the rc4 inventory partial). The head units
   enter at their planned widths (`validation.components_drawn.head_unit_widths_um`), and the head's gap mouths (the
   gaps of the rows below NM) are not in the inventory, so both are held only by the DRC at the floor F_min (item 31;
   the check's own text and `auditor.drawn_widths` still call every width measured). No shipped outcome moved,
   because on every shipped design the narrowest channel is either one that rc4 already measured or the gap.
   Strict mode refuses (exit 5) or fails (exit 2) a design above the limit, so a strict PASS never carries a channel
   narrower than w_min. Evaluate and reproduce modes draw such a design and report it (kind `fabrication`; the verdict's
   ending "ABOVE THE PROCESS ASPECT-RATIO LIMIT ... reported, not refused or failed in reproduce mode";
   `assessment.over_aspect_ratio`, `api.Result.over_aspect_ratio`, and `status.aspect_ratio_of_CIF_channels`). G11's
   fold itself is not implemented, so where spec v5 would fold a narrow side channel, opencif refuses the design
   (strict) or reports it. The five fine-floor examples (`m13`, `m15`, `m16`, `m17`, and `gifford_physics_linear`),
   whose 4.3 to 4.5 µm floor at 80 to 100 µm depth already assumes a finer process (15 to 23 : 1), declare it as
   `process.aspect_ratio_max: 25`, which is a stated process capability rather than a fact about any foundry. The two
   Gifford design-level reproductions (`gifford2014_fig4_design` at 16.5 and `gifford2014_fig6_design` at 32) keep their
   geometry in reproduce mode and are reported, whereas their as-fabricated presets use the 8 to 20 µm floors (10 : 1 or
   below). The inlet stem and filter are checked against the floor (DRC) rather than against AR_max. The limit is
   applied to the fluid channels (the ridges of an SU-8 mold). The walls between channels (the mold's trenches) have a
   limit of their own, which opencif checks only through F_min, by the DRC (a review of rc4). The future work is G11's
   fold geometry and a wall aspect-ratio rule if a process states one.

## B. The 1-D network model

8. The junction terms are 2-D fits (κ_j for bends, c_T and its fillet term for head T-junctions, and, since rc5, the
   outlet stretch-end term c_exp of 15a), placed between their deep and shallow limits by a depth proxy. The weight is
   the curved-duct solution's (`schedule.depth_lambda`, which gives where the exact conductance of a tight reference
   bend, r_m = w, at the channel's depth and width falls between its planar and Hele-Shaw values), and it is carried
   over to these junctions rather than solved for them. κ_j and c_T matter most at a head's seam rows, and c_exp matters
   most at the last rows. The record gives the model band (the rows that the terms' range moves, with their share of F),
   in which the depth weight is moved by ±max(weight, 0.1), the T-junction term is widened by ±(a quarter of its fillet
   part + 0.02 w), and c_exp is scaled by ±30% (rc5). The band is a sensitivity scenario of a stated size. It is neither
   a calibrated interval nor a proven bound.

   After the rc5 reviews, its c_exp part was checked in 2-D against three things that it stands in for
   (`provenance/rc5_scripts/`). The first is (i) a depth proxy that is wholly wrong. On all 21 designs with the term,
   c's deep and shallow limits both lie inside the band (they span −31.9% to +54.9% of the recorded c on Xia's centered
   rings and −14.9% to +1.3% on Dinh's flare, and the narrowest margin is 1.3 points, on the three rings at width ratio
   2.19; `band_coverage.out`). This does not hold in general, because one limit falls outside the band for a one-sided
   step at width ratios below about 2.25 (at small weights, and below 1.07 at weight 0.46) or above about 5, and for a
   centered step below about 2.9 (at small weights) or above about 7.7. The second is (ii) the drawn flare's 25 µm
   tangent arcs, where the tables have sharp corners. Their effect is −0.7% (deep) and −0.2% (shallow) on Dinh V1, the
   same on two meshes (`rounded_flare.out`). The third is (iii) the teardrop lead's first quarter arc (centerline radius
   2 lead widths), which starts at the stretch's end, where the tables put a straight wide channel. It adds +8.0%,
   +6.4%, and +4.7% to a one-sided step's excess in the deep limit and +26.2%, +22.6%, and +18.6% in the shallow limit
   at width ratios 2.4, 2.9, and 3.87, the same on two meshes (`step_bend.out`), so counting the arc would raise c by
   8.0% to 14.4% on the 9 teardrop steps. That rise lies inside the band's upper side (+22.8% to +32.3%), although it is
   a one-sided bias rather than a symmetric uncertainty. The ring leads' turns (with their own radius and angle) and the
   arc after Dinh's flare have not been solved. Release rc5 also corrects rc4's band, which never moved the outlet
   leads' arc terms, because it read their arcs from the measured outlet record, which does not carry them. On the
   shipped designs, that term is small beside the others. Where an outlet does not change width, it moves the last 1 to
   6 rows of the seven Gifford presets with separate outlets by at most 2.1% of D, and the last 6 rows of the 19 µm
   Gifford-targets example (whose record gives the band's largest change, 14.4% of D, at a seam row rather than at the
   last rows). The 2-D tables end at width ratio 8 (steps) and at 4 and 45° (flares). Beyond them, a step's deep value
   is held at the table's end (its shallow value is the closed form), and a flare's deep and shallow values are both
   held there. Flares between 1° and 5° interpolate toward the 0° column (no excess). No shipped design reaches a
   table's end, because the largest ratio is 3.87. A user design beyond one is held at it without a message (rc6: a
   named status, and in strict mode a refusal, outside the solved range; a review of rc5). Neither appendix discusses
   these terms. The future work is a 3-D Stokes calibration of κ_j, c_T, and c_exp at real aspect ratios, with the
   lead's first arc (`tools/cfd_kit.py` writes the seam and U-turn windows and, as yet, no outlet window).

   Since rc7 (DINH-1), a square quarter turn (`meander.turns: square`) carries the corner term c_corner, the excess over
   the centerline through the corner point, from the same two 2-D limits (`provenance/rc7_scripts/corner/`, with
   Hele-Shaw finite volumes and planar Stokes by the dissipation integral at 40 and 80 cells per w, first-order
   extrapolated; the solver reproduces rc6's sharp-corner values and, on the arc, rc6's exact curved-duct solution to
   0.004 w). The term is placed by the same depth weight and banded by ±(¼ |c − c_sharp| + 0.02) w per corner, where
   |c − c_sharp| is the inner rounding's effect. Its table runs over the inner rounding r_in/w from 0 to 0.44 at r_out/w
   0.25 (the outer rounding moves c by under 2 × 10⁻⁴ over 0 to 0.31), and beyond r_in/w 0.44, c is held at the table's
   end (`model_scenarios.bend_model.square_corner.held_at_table_end`). On Dinh V1 (r_in/w 0.36), the band moves 121 rows
   beyond ±1% (111 head rows, 103 to 213 or 104 to 214 at the band's two ends, and the last 10; 4.3% of F), against 11
   rows with rc6's arcs (1 head row and the last 10). The corner is therefore the head's least certain element. The
   drawn-mask measurement takes the roundings' area off before it converts a unit's area into length, because the corner
   term counts that area. The term c_corner is the glossary's L_eq (the "bend-equivalent length" that the glossary lists
   among its External Inputs as supplied by this work; spec v5, section 14.6), which spec v5 asks to be
   "calibrated once at low Re". Here it is a 2-D computation in the two limits with the depth weight rather than a 3-D
   calibration, which remains future work with κ_j, c_T, and c_exp. Its sensitivity to the inner rounding is about −0.9
   per unit r_in/w (−0.032 per µm at w = 28 µm, so ±2 µm of r_in moves a corner by ±1.8 µm of length), and the band's
   ±0.10 w per corner covers about ±3 µm. Relative to the centerline through the corner point, the term is negative
   (−0.24 to −0.56 w), so spec v5's "typically 2–6 %" excess and "under-resistant" folds (S:407, S:757), which are
   statements about U-turn folds measured along their own path, do not transfer (proposed spec erratum D105). Spec v5's
   G4 (a constant cross-section along every traversal) does not hold at a square corner, which is why the term exists,
   whereas G12 (turns counted from the drawn polyline) holds by construction (`geometry.unit_path` records each corner).
   Both radii (7 and 10 µm) exceed the process's corner limit R_min, which spec v5's G7 requires of inner corners and
   which here is the generator's default min(4.5, 0.25 × floor) = 4.5 µm (G; the Dinh presets set none). The head
   smoothing rounds every corner at the 5 µm fillet at least (code review R1). Square quarter turns lie outside spec
   v5's own turn geometry (semicircular U-turns).
9. The gap entrance coefficient c_ent = 1 (external) scales every gap resistance and so every f. The 3-D gap-mouth
   effect on the per-gap f_eff, which App. B C4 CQ1 calls a "sibling" of the yield-side residual (an uncertainty on √c
   on the scale of about 5%), is not modeled. The record carries the scenario "mouth R_gap ±20 %" (head mouths and the
   rectilinear rows only; the linear gaps' c_ent is not varied). The future work is a 3-D calibration.
10. The side channel's open inner wall. The 1-D rule treats a straight side channel as a closed duct. A Hele-Shaw
    estimate puts the side resistance 1.5% to 7% low (w_s from 100 to 20 µm), and the error varies along the device, so
    κ_g cannot absorb it. This limitation is framework-level, because the published designs share it. The future work is
    a 3-D solve of a few rows.
11. Dean (inertial) terms at bends are not modeled. Dean's circular-pipe leading term, 0.03058 (K/576)² with K = 2 Re²
    a/R, is about 1.5 × 10⁻⁴ of a bend's resistance at Re = 10 and a/R = 0.2, and 3.7 × 10⁻³ at a/R = 1, which no
    drawable arc exceeds (here a is the half-width in the plane of the bend, and R ≥ a is its centerline radius). A
    channel of depth h and width w carrying Q_x has Re = 2ρQ_x/(μ(w + h)). Because every channel carries at most the
    inlet flow and is at least the DRC's tested width (F_min − 0.1 µm) wide, Re < 2ρQ/(μ(h + F_min − 0.1 µm)). The
    record gives that bound and, where it is ≤ 10, the Dean term at it with a/R = 1 (`Re_any_channel_bound` and
    `dean_leading_term_bound`). On 38 of the 50 drawn presets and examples (the three Fig. 4 trace variants excluded), the bound is ≤ 10 and the Dean term ≤ 3.5 × 10⁻³ (the
    refused `wbc_xia_process_refused` is not drawn, so its DRC-conditional bound does not apply). On the other 12
    (Gifford Fig. 4 at 30 µL/min and 80 µm deep and Fig. 6 at 55 µL/min, their aligned readings included, the 19 µm
    Gifford-targets example, and the four Xia RS-7 presets), the bound is 10.0 to 11.9 (Xia 171), whereas their gated
    center-channel Re is 5.6 to 6.7 (Xia 61, outside the envelope). The recorded term is an order-of-magnitude screening
    estimate rather than a validated bound on a rectangular duct's bend error, because the coefficient is a circular
    pipe's, the per-bend a/R is not computed, and no 3-D solve stands behind it. With a taken as half the hydraulic
    diameter, D_h/2 (which approaches w in a deep narrow channel), the term is up to four times larger.
12. Some parts are taken from the design rather than measured on the drawing. They are the radii of head, U-turn, and
    outlet arcs (their widths are measured), the head rows' center-path lengths and the lengths of the generator's
    straight segments, and off-chip tubing. The record's notes list them. Since 1.0.0, the physics audit measures the
    straight rows' lengths and the post length across the flow (ℓ_g, which every gap resistance takes from the design)
    on the photomask file itself. Its own ruler reads the row pitch and the wall thickness from the posts and holds the
    generator's measurement to them (R0), and R6b requires every straight-row post's extent across the flow to be ℓ_g
    within 0.5 µm. The head units, the U-turns, and the outlet are not re-measured by the auditor. The generator's own
    measurement of a head unit takes its length as its area over the median of width sections at 25, 50 and 75% of its
    straight segments longer than three widths, so a local narrowing changes the measured resistance little and with
    the wrong sign (a 10 µm stretch narrowed to half width, seeded on `dc_D5_F0p3`, measured −0.19% against +2.4%), and
    a unit without such a segment takes its planned width (48 of the 219 head units of `dinh2024_v1_depicted`). A
    smooth width defect of a head unit above the floor is therefore seen only weakly by Q1 (item 31).

## C. Sensitivity and fabrication

13. Seam and turn amplification. A meander head's seam rows and the rows next to a tight U-turn amplify element errors
    about 70-fold, so a few 10⁻⁴ of resistance move their cutoff by percent. The record carries scenarios on the design
    network (`validation.robustness_design_network`), which are a head-unit junction length of ±w/2, a uniform width
    bias of ±0.5 µm, a mouth R_gap of ±20%, and one seeded draw of 10 nm per-segment side-width noise (a single
    realization rather than a distribution). There is no depth scenario. The perturbation sizes are illustrative
    assumptions rather than a fabrication process's measured tolerances. Their size on the shipped designs is given in
    the RESULTS.md columns, with counts over both signs. A ±0.5 µm bias moves 17 to 223 gaps of a device without a head
    outside D ± 1% (up to 21% of N, with cutoffs of 0.84 to 1.16 D, worst at a sealed start's first rows) and 82 to 404
    gaps of a headed device. For a headed device, the counts include its declared rows, which are outside at nominal,
    and they reach up to 23% of N, on `dc_D3_F0p3` and `dc_D4_F0p3`, whereas the 24% of `dinh2024_v1_figure_aligned` is
    mostly its declared head rows, the bias itself moving about 4%. In a headed device, the bias is worst at the seam
    rows, where under −0.5 µm one row of 24 of the 33 headed designs loses its cutoff (reverse flow or f ≥ ½), and the
    head-unit junction length −w/2 does the same to 1 to 3 rows of 29 of them. The noise draw moves 0 to 2,005 gaps. The
    statement "D ± 1 % at every gap" therefore holds for the nominal geometry rather than for a fabricated chip. The
    future work is a fabrication tolerance analysis, including a tolerance budget per record (the bias at which the
    share of F extracted within ±2% or ±5% of D falls below 95%).
14. Per-gap cutoffs are far more sensitive to sub-micron side-width errors than the device-level separation is, so
    report the share of filtrate extracted off-cutoff (the record's `share_of_F_within_1%`) in addition to the worst
    row.
15. Off-chip tubing enters the outlet closure as a design value. The record carries tubing scenarios and dD_last/d ln
    R_ext.

15a. The last rows rest on the outlet balance (rc4; rc5 models the stretch end). The last rows' gap flows are a small
    difference of the two outlet pressures (on Dinh V1, the last gap drops about 0.06 Pa, whereas each outlet sits about
    110 Pa above the sink), so their cutoffs move with any error in the side outlet's resistance R_s,out. Over the 50
    drawn designs, d(D_N/D)/d ln R_s,out is −26 to −175 (the record's `dDlast_dlnRside_out`). The sensitivity belongs to
    the closure rather than to opencif, because any CIF device whose last gaps are fed by the outlet balance has it, the
    papers' devices included. This is one reason why the measured F (Dinh's 0.505 against the outlet identity's 0.559)
    is the only outlet-level check that the literature offers. Three parts of R_s,out are outside the 1-D duct rule:
    - The first is the width change at the end of the outlet stretch, where the side channel meets a lead of another
      width (every teardrop and ring outlet, and a separate outlet only where w_s(N) < 40 µm, which no shipped design
      has). Release rc4 did not count it, and its documentation sized it wrongly as "a fraction of w", illustrated by
      0.05 to 0.25 w. The 2-D solves of rc5 (`provenance/rc5_scripts/expfem.py`) put it, on rc4's masks, at 0.42 to 0.59
      of the narrow width for the one-sided steps of 18 shipped designs and at 0.14 for the centered narrowing of the
      three Xia RS-7 rings, which is 0.34% to 3.49% of R_s,out (`provenance/rc5_scripts/rc4_step_excess.out`). On Dinh
      V1, it is 0.513 w, or 3.49%. Added to rc4's Dinh network, which did not count it, it puts the last 20 rows outside
      D ± 1%, 8 of them without a cutoff (reverse flow), and moves F from 0.5593 to 0.5521, beyond Q2's 10⁻³
      (`provenance/rc5_scripts/rc4_step_effect.py`). Since rc5, the closure counts it (`schedule.outlet_exp_c`; item 8)
      with a deep (planar Stokes) and a shallow (Hele-Shaw) value for the geometry (one-sided step, centered step, or
      flare; width ratio; flare angle), placed by the depth weight, so the lead's drawn resistance is lower by that
      excess and the term's range is in the model band.
    - The second is the 3-D flow at the outlet bifurcation (the glossary's measured η = 0.903 is a 10% effect on F),
      which is not modeled.
    - The third is the port, which is not modeled.

    Every record reports both side outlets at ±1% of their resistance, an illustrative size like the other scenarios
    (the record's `model_scenarios` `outlet_side_resistance_±1%`, and the model band check's "outlet balance"). Over the
    50 drawn presets and examples (the three Fig. 4 trace variants excluded), ±1% moves 12 to 59 rows outside D ± 1%, always the last rows (within the last 59, and Xia RS-7's
    within its last 80), by up to 0.30 to 1.06 D (0.02 D on Xia RS-7, whose off-chip tubing dominates its outlet),
    carrying 0.07% to 2.9% of F. The stretch-end term's own range (the depth weight moved, and the term scaled by ±30%,
    a sensitivity scenario whose span item 8 describes) moves the last 6 to 58 rows of the 24 designs whose outlet
    changes width beyond ±1%. These rows appear in the record's junction band, which moves every junction term at once
    and whose largest cutoff change over all its rows runs from 0.004 D on Xia RS-7 to 0.76 D on
    `gifford2014_fig6_aligned`, and on the three Fig. 6 presets, 10 of those rows lose their cutoff at an end of the
    range. The band's share of F, seam rows included, is 0.02% to 1.9% on those designs (4.3% on `dinh2024_v1_depicted`,
    whose square-corner term adds 111 head rows, item 8). The seam rows' band is unchanged from rc4. The record's
    junction band names the last rows beside the seam rows, and since code 56e1980d8a064b35 the record lists every row
    of a scenario, whereas rc4's listed 40. Thus, "D ± 1 % at every gap" holds for the nominal 1-D outlet rather than
    for the last rows of a real device, whereas the device-level separation is barely touched. The future work is a 3-D
    solve of the outlet (calibrating c_exp and the bifurcation) or a terminal design that makes the last rows
    insensitive.

15b. The outlet leads' shape. A teardrop, separate, or ring outlet runs the last side channel straight for `L_out` and
    then changes its width into the lead. By default, the change is a right angle filleted at `fillet_um`, with the
    inner wall continuing and the outer wall stepping out (a one-sided step; on rc4's Dinh V1, 44 → 119 µm at x
    138.68 mm), or, for a lead narrower than the side channel (Xia's ring), a centered step. Release rc5 adds a flared
    lead, which is opt-in (`outlet.flare_deg`, `outlet.flare_start_um`) and available on teardrop and separate outlets.
    The flare is refused with the ring, with `outlet.stretch_um`, which the flare replaces, and with
    `physics.fail_safe_margin`, whose re-solve rebuilds the march and so moves the last side channel's width on which
    the flare was planned (since rc6, the refusal's message says so). Since rc6, it is also refused with a
    reproduce-mode `reproduction: {L_out: …}`, which rc5 ignored with a flare on a linear layout. With a flare, the side
    channel runs straight from the last gap's center to the flare's virtual corner at `flare_start_um` (or, where the
    last leg ends beside a loop, to the longer of that and the loop's clearance), and then its outer wall diverges at
    `flare_deg` per side, with 25 µm tangent arcs at both kinks, to the lead width, whereas the inner wall stays
    straight. A lead no wider than the side channel (a separate outlet with w_s(N) ≥ 40 µm) draws no flare and keeps its
    stretch. The record's `design.outlet.L_out_um` is the drawn stretch. The widening is integrated by the exact wall
    rule (the drawing and the integral use the same polyline), and its excess over that integral is the stretch-end term
    of 15a (0.08 w at 14.5° and width ratio 2.4, against 0.49 w for a step of the same ratio). The two Dinh V1 presets
    draw it with the values measured on Fig. 2B's fourth panel, 14.5° per side after a straight 65.4 µm (rc6, at the
    50 µm bar read with its anti-aliased ends, whereas rc5 drew 67 µm; PRESET_PROVENANCE; `tests/test_figure_anchors.py`
    checks the drawn walls), so they now reproduce the side channels' outlet end (REPRODUCTIONS). The lead width that
    the flare reaches (105.75 µm) is the generator's (A), because the panel does not show it. The angle is measured from
    each wall's own straight direction, and that direction is the least certain input. Fitted over six row ranges of the
    native raster, the straight walls tilt 0.28° to 0.40° (left) and 0.51° to 0.73° (right), giving 14.43° to 14.55° and
    14.38° to 14.60° per side (a review of rc5 found the right wall at 0.09° to 0.35°, or 14.7° to 15.0°, on an image at
    four times the native resolution, which was not reproduced on the native raster). The preset's 14.5° lies inside
    both native-raster ranges (`provenance/rc5_scripts/measure_dinh_centre.py`). The panel also shows the center
    channel, 49.7 ± 0.2 µm outline to outline after the last gap (52.8 µm between the post columns, printed 52 µm),
    which begins to widen 191 µm after the last gap's center at 0.25 µm per µm, to 58.5 to 58.9 µm at the panel's lower
    edge, 228 µm after it (the same script's picks at the corrected bar, rc6, whereas rc5 read 50.7, 53.8, 195, 59.6 to
    60.1, and 232 µm at 46 px; every pick is in `dinh_centre_measurement.json`). opencif keeps the center stretch at
    w_c(N) to its end, and the center lead is its own (A), because a widening seen over the panel's last 37 µm does not
    give the rule that the lead follows beyond it, and the outlet closure would absorb it in the center lead's length.
    Every other outlet keeps its step, now with the step's term in its closure.

## D. Generator coverage

16. The tapers are constant and constant-shear only, with no Ψ (canonical) taper of App. A R13(b) and no user profile.
    Since rc7 (XIA-1), the constant-shear taper holds either the rectangular duct's maximum wall shear
    (`shear_formula: exact`, CIF_v2's C(α), spec v5 S16-b) or the parallel-plate shear 6Q/(a b²) (`parallel_plates`, w_c
    ∝ Q_c). The second is a departure from spec v5 ("NEVER 6") that is allowed in reproduce mode only. Which one a
    paper's device held is a reading of its figures (for Xia, register XIA-1). Xia's text cites
    "formulas for rectangular channels based on Yang et al.", whose Eqn 3 the glossary identifies as the resistance
    series (G:2036 to 2038), and does not state the shear coefficient. Every record reports the peak wall shear via C(α)
    along the center channel (`validation.gates.shear_peak`, as spec v5 S:606 asks). Under the parallel-plate taper of
    `xia2016_rs7_aligned`, it rises 17.4% along the taper and is 1.34 to 1.57 times the parallel-plate value. (Spec v5's
    constant-shear ODE, S:456, holds only for w_c <= h, and the code solves the shear numerically and does not use it.)
17. Composition is limited to parallel arrays with per-element outlets. Merged outlets, series stages, and cascades are
    not implemented (merging outlets changes each element's outlet closure).
18. Dinh's filter is drawn in his absolute coordinates and fits only his three-arm head. The ring outlet is parametric
    (linear and paperclip layouts only). Its center lead is fixed by the ring, so at a high F it can be unable to close
    the design, because its side lead would need less resistance than the side channel's straight run already has, or a
    lead too wide for the ring's turn. Since rc4, such a design is found before drawing, by planning the outlet on every
    leg count that the layout would try, and it is refused by name (exit 5, "REFUSED (outlet)", in `check` and
    `generate(check_only=True)` alike), whereas rc3 failed inside the drawing code with a root finder's text. The
    pre-check plans the outlet at the head width that the spec gives (or the default). When the fitter widens a packed
    head afterward (to seal gap 1 or to fit), the outlet is planned again while drawing, so an outlet that closes at the
    nominal head width and not at a widened one is refused by the build (exit 5, the same refusal record) rather than by
    `check` (a review of rc4; the 22-case build lattice met no such case). The flared lead (15b) is not drawn with the
    ring, and no shipped spec or test builds it on a separate outlet or on a layout whose last leg ends beside a loop
    (only the two linear Dinh presets draw it). The plan-stage lattice (`tests/test_lattice.py`) has 42 plans (4 named
    processes × 3 outlets × F 0.1, 0.5, and 0.9 on a 25 × 75 mm paperclip, plus other layouts). In it, a ring at F 0.9
    is refused on the generic, Gifford, and Xia paperclips and on a linear generic device, and it plans on Dinh's
    process. Every ring plans at F 0.1 and 0.5 where the layout itself can be laid out. At F 0.1, the paperclips of the
    generic and Gifford processes are too short to fold, and Dinh's process reaches F inside its head, and both are
    named spec errors before drawing (rc4), whereas rc3 raised the second as an IndexError. The shipped rings are drawn
    at F 0.5 and on the papers' measured rings. The suite plans the lattice without drawing it, and a full build of its
    22 cases was run once for rc4's verification (`provenance/rc4_scripts/build_lattice.py`). The future work is a
    pairwise build lattice that is run routinely.
19. Sealed heads at small D are long. A sealed head at Dinh's targets draws a device 167 × 4.2 mm (a 175 mm chip),
    against 148 to 152 mm for Dinh's own unsealed head as reproduced here, and the papers' own heads are drawn unsealed
    only with gap 1 declared. Hairpin heads cannot seal gap 1 at most targets.
20. Head capacity planning uses geometric lengths (rc6), less a square turn's corner term c_corner w (rc7, code review
    R7, because without it a square-turn unit's capacity was 25.5 µm over its effective length at Dinh's full reach).
    The U-turn's arc term and the T-junction's term are still left out, so a unit at full reach can draw a few µm longer
    than planned. At Dinh's head, the capacity exceeds the effective length by about 9 µm (8.8 µm with arc turns and 9.7
    with square turns; the review's f6_capacity.py), which is well inside the 15 µm lateral margin, and the envelope
    check (outer edge against W/2 + 0.5 µm) catches more. No shipped spec reaches its capacity (Dinh's largest taper
    unit needs 2,977 µm of 3,002).
21. A build takes 1 to 30 min per device on rc4's machine (2.5 to 3 h for all 47 specs on two cores, rc4; here up to
    53 min alone and 85 min beside the suites). The first full rebuild of rc5 took 4 h 22 min, and its final one took
    5 h 30 min (alongside seven test suites), on two cores on a machine running about twice as slowly (rc4's code, timed
    again on it, ran m15 in 90 s against its recorded 46 s). After the rc5 reviews, the two Dinh presets were rebuilt
    (their notes only, four times) in about 12 to 13 and 17 to 23 min each, one at a time (the last beside a test-suite
    lane).
22. The realized cutoff is a model prediction on the drawn geometry rather than a measured particle cutoff, and no
    device in this repository has been fabricated from these masks.

23. User processes below the tested gap and floor (rc6, from Reviewer A's E5). The named processes (`api.PROCESSES`)
    have gaps of 16 to 20.5 µm and floors of 16 to 19 µm, and the sweeps, the matrix, and the plan-stage lattice stay
    within that range (the Gifford presets' lower floors, 4.5 to 15 µm, draw no head). A process of 12 µm gaps and a
    10 µm floor (60 µm deep) with the narrow head drew its head without gap openings in rc5, and the drawn-mask
    measurement crashed (`ZeroDivisionError` in `measure2`, now `opencif.checks.measure`). The cause was the head and
    outlet smoothing (a closing and then an opening of radius `process.fillet_um`), which erases every channel no wider
    than twice the radius, as the 5 µm default did at a 10 µm floor. Since rc6, the default is min(5, 0.4 × the
    narrowest drawn feature, 0.3 × the narrowest head unit) (3 µm there), an explicit fillet at or above half the
    narrowest feature is refused before drawing, and a feature that the measurement cannot find on the drawing is a
    named FAIL ("every planned feature found on the drawing") and never a crash (`tests/test_negative.py` 13). In such a
    FAIL, the record's measured widths of the missing features are the plan's design widths, put in so that the
    measurement finishes (`validation.features_found` lists the features, and `by_kind` counts them). The fields built
    on them (`fidelity.gap_width_dev_max_um`, the loops' and leads' narrowest widths, `drawn_parameters.G_meas_min`, and
    the aspect-ratio inventory) are therefore not measurements for those rows. At the default, the same device draws
    every feature and builds to the end. It fails Q1 and Q2 only on three rows at the head-array seam (rows 221 to 223
    of 1,512; cutoffs of 0.937 to 0.989 D; 0.18% of F), which is the seam amplification of item 13 at a process below
    the tested range, and the failure is reported by name (exit 2; `provenance/rc6_scripts/e5_default_check.py`;
    `tests/test_negative.py` 13 checks that every feature is found and that the topology is as planned). Below the
    tested range, then, the generator either draws the device or says by name what it could not draw, and we do not
    claim that every such process yields a passing design.
24. What the micrograph traces show, and what they do not (rc6; `provenance/rc6_traces/`). Every micrograph panel of the
    three papers that shows post rows, loops, or turns, and Dinh's head, is traced (Xia's head, Fig. 3b, is measured but
    not registered) with the shipped rc6 mask at the panel's own scale. That scale is a bar read with its anti-aliased
    ends or, for Xia's Fig. 3 panels b to g, whose bar is panel i's, their printed widths, with panel f at its own
    period. Each panel is registered twice (features and chamfer) and judged against its nulls and a scale scan. The
    traces show the drawn geometry where the papers' devices are visible, on Dinh's four panels at medians of 0.44 to
    1.29 µm (the outlet's widths within 0.5 µm of the print), on Gifford's Fig. 4 inset at 1.21 µm and Fig. 6's four
    legs at 1.4 to 1.9 µm each (and on one map with their turns at 1.47 µm), on Xia's panels at 2.3 to 3.8 µm (relief
    contrast), and on Xia's whole device on Fig. 3a. These are rc6's medians, on rc6's masks and window; the shipped
    masks' medians (rc7 and RC8, under the window rule below) are tabulated in `provenance/rc7_traces/README.md`
    (Dinh 0.44 to 1.40 µm; Gifford Fig. 4 leg 4 0.61 µm nested, 1.10 µm design; Fig. 6 legs 0.51 to 0.64 µm aligned,
    1.09 to 1.87 µm design; Xia's claimed panels 2.02 to 2.78 µm). The pitch and post sizes of the presets were measured on these same
    figures, so their agreement is consistency. The same holds for the layouts (Gifford's leg spacing and leg axes, and
    Xia's Fig. 3a leg axes), which were read on the figures that they are traced against, except for where along the
    flow Fig. 6's turns start (the mask's 5.7 µm downstream of the image's). It also holds for any width that sets a
    panel's scale or places it along the device (Xia's printed 300 / 150 µm and its loop curvature; the post columns
    that place Dinh's panel 3) or that was read on the panel (Dinh's flare). The tests are the legs' ends and identity
    and the widths that the schedule sets where nothing read on the same panel fixes them. These are Dinh's side
    channels on panel 3 and its outlet widths against the print on panel 4, Xia's taper on panels d and e and its last
    side channel on g (since RC8, panel d on both masks and g on the depicted one are recorded rather than claimed; A7,
    D5, register XIA-9), and the side channels at Gifford Fig. 6's leg ends. The limits in `trace_anchors.json` were set
    with the results in view (regression limits). What the traces do not show falls into six cases, (a) to (f).

    (a) Dinh's transition blocks read about 2.2 µm shorter at each end along the flow than drawn (V3 block_along −2.10 /
    −2.47 / −2.23 µm, center / left / right). The cause is either the outline bias of a relief image or blocks drawn 2
    (s − GS) long where the device's are shorter, and the two cannot be separated on one panel.

    (b) Gifford's Fig. 4 neighboring leg lies 818 µm from leg 4 on the inset (the generator's equal-leg rule gave 969,
    and Gifford prints no spacing). The presets of rc6 take the inset's spacing, so leg 2's posts lie on the image's by
    construction, whereas its side channels there read 17 to 18 µm wider per side than the rows that the mask places at
    that point (rc6: which rows pass the inset, the inset did not fix; RC8: on the nested reading they are leg 2's rows
    1,730 to 1,758, item 28). The claim's 90th percentile (18.2 µm) was the metric's search window rather than geometry
    (rc7, register GIF4-2), and with a half-window of at most half the gap, which is the rc7 window rule, it is 2.5 µm
    on the design mask. The aligned reading of rc7, reading B (legs in path order; since RC8 the rejected variant
    `gifford2014_fig4_pathorder`), put the neighbor, read as leg 3, on the image's walls through a fitted leg 4 end, and
    RC8 rejects it by the taper direction (V13; register GIF4-1, item 28).

    (c) For Xia's center-channel taper (rc7: the parallel-plate reading of XIA-1 and item 16), the panels read 300 µm at
    leg 1's first posts, 220 µm where leg 3 leaves its loop, 186 µm half way down leg 3, and 150 µm on legs 5 and 7,
    against the assumed constant-shear taper's 300, 233, 205, and 150 µm. According to PRESET_PROVENANCE, Xia prints the
    ends and the principle, a relatively constant maximum shear, rather than the profile, whose form is an assumption
    (A). Thus, the device narrows faster than the assumption over legs 2 to 4.

    (d) For Xia's length, at the pitch that the panels show, Fig. 3a's legs hold 1,190 rows where the schedule places
    1,079, so the whole-device registration on Fig. 3a puts the drawing's far loops 5.7 mm and its outlet ring 5.4 mm
    beyond the mask's (`xia_fig3a_device.py`). The variant with D 6.5315 µm (× 0.9515, F unchanged) fills the legs (its
    far loops 0.09 mm and its outlet ring 0.19 mm from the drawing's, and panel g 0.65 mm from its box, 0.46 mm of which
    is the panel's own offset from its box). Under the model's premises, this is the cutoff that the drawn device
    implies, a conditional like Dinh's length. The variant with F 0.94072 (panel g's 624.6 µm last side channel) and D
    6.5705 µm (Fig. 3a's legs at that F) matches both (panel g V5 across 1.005; far loops +0.09 mm, ring −0.19 mm).
    However, panel f, a quarter of the way down leg 7 (recorded at its own 1.774 µm/px, not claimed), shows side
    channels 29 to 33 µm per side narrower than that variant draws (20 to 24 µm narrower than the preset, and about 7 µm
    less at 1.798 µm/px). If panel f's scale holds, the device's side channels widen faster along the last leg than the
    model's, so no F or D matches both panels (the schedule's form, as in c). In rc7 (register INV-4, open), on the
    shipped `xia2016_rs7_aligned` (parallel-plate taper, D 6.4159 µm), panel f's walls lie 39.9 / 36.0 µm inside the
    mask's (left / right, at panel f's 1.774 µm/px, and about 7 µm less at 1.798), which is more than on rc6's variant,
    whereas panel g matches (−0.3 µm) because F is fitted to it. Neither the taper option nor a parallel-plate
    side-channel law (which would need F 0.893) explains the difference.

    (e) On Gifford's Fig. 6 (new in rc6), the side channels at the ends of legs 7, 5, 3, and 1 read 5.4 to 10.0 µm per
    wall wider than drawn (145.5 / 114.3 / 87.3 / 58.9 µm against 139.9 / 108.7 / 80.7 / 49.3; per leg 5.5 / 5.6 / 6.6 /
    9.6 µm, growing as the channel narrows). The assumed gap (Fig. 4's 20.5 µm, whereas this image reads 17.1 along the
    flow), the outline bias of straight walls against rounded tips, or the device could each cause the difference. The
    generator's uniform layout (723 µm, turns of radius 362 µm) did not match the image (legs 588 µm apart; turns about
    one center 1,745.7 µm beyond leg 1, its SD 2.2 µm in the common fit and 26 µm between the walls fitted one by one,
    radii 1.68 to 3.17 mm), so the presets of rc6 take the image's leg axes. The axes of legs 2, 4, and 6 come from
    assuming semicircular turns, because only 8° to 23° of each turn is in the frame. The turn walls then lie a median
    1.4 µm from the mask's edges, which is consistency because the axes were read on those walls, and the turns' start
    along the flow, which the walls do not set, lies 5.7 µm downstream of the image's. On the image, the last leg's
    array ends 5 rows before the others' (item 25).

    (f) For Dinh's head (Fig. 2B, first panel) on rc6's masks, the generator draws every meander unit's turns as arcs of
    centerline radius (channel + wall)/2 = 27.5 µm (outer wall 41.5 µm, inner 13.5 µm; `opencif/geometry/shapes.py`,
    `unit_path`), where the micrograph shows square outer corners. At a corner, the two depart by up to about 17 µm,
    beyond the metric's ± 8.5 µm window near the apex, so those stretches drop out of the panel's median and 90th
    percentile (gray on the residual map), whereas the flanks, at −4 to −6 µm, are counted. The head's resistance is
    computed on the arcs drawn (with the bend correction), so the model is consistent with the mask, and the shape is an
    assumption (A) rather than a reading. Release rc7 draws the depicted head's turns square (`meander.turns: square`;
    register DINH-1; item 8). Arcs stay the default, and `dinh2024_v1_specified` keeps them (the text prints no corner
    shape).

    The rc7 traces (`provenance/rc7_traces/`, the eight shipped masks that the figures show) keep rc6's limits, with
    exceptions that were all set with the results in view. The exceptions are Dinh panel 1's V5 across (0.975, 1.015,
    where the square head prefers 1.011 and rc6's arcs 0.988) and its null test (chamfer, the metric null ratio being
    1.17), the Fig. 4 design's p90 (3.5 µm), the recorded Xia f's V5 across range on the parallel-plate mask (0.94,
    around its 0.950), and the new masks' claims (Xia e on that mask promoted to a claim with a chamfer null ratio of at
    least 1.15, and its panel d median limit of 3.5 µm). In rc7, 21 null ratios lie below the package's 1.3 criterion,
    14 of them on 11 claimed panels. Six claimed panels pass V2 on a chamfer ratio of 1.05 to 1.28 (Fig. 6 design legs
    1, 3, 5, and 7 at 1.17 to 1.28; Xia g on the depicted mask at 1.048; Xia e on the aligned mask at 1.21), and the
    metric ratio is below 1 on Xia c (0.81, both masks) and on the aligned mask's Xia e (0.53), where a
    half-pitch-shifted mask scores better than the registered one. The width readings do not depend on the nulls,
    although the residuals quoted for those panels are weak evidence (`trace_anchors.json` null_flags). Since RC8 (A7,
    D5), every record states its coverage and its censored null ratios. Xia Fig. 3d (both masks) and 3g (the depicted
    mask) are recorded rather than claimed (coverage 0.42 to 0.48, censored null ratios 0.94 to 1.00; register XIA-9),
    so panel d has no median limit and Xia g on the depicted mask is not a claim. With the censored ratios counted, 33
    null ratios lie below 1.3, 16 of them on 9 claimed panels. Five claimed panels pass V2 on a chamfer ratio (Fig. 6
    design legs 1, 3, 5, and 7 at 1.17 to 1.28; Xia e on the aligned mask at 1.21), and the metric ratio is below 1 on
    Xia c (0.81, both masks) and on the aligned mask's Xia e (0.53).
25. A paperclip's last leg takes the rows left over (found on Gifford's Fig. 6, rc6). The layout divides the path into
    equal legs and ends every middle leg at a whole row, so the last leg's array can end several rows beyond the others'
    ends. On Fig. 6 at the generator's defaults, the middle legs hold 869 rows and the last 875, so its array ends 6
    rows, 308 µm, beyond the other legs', whereas Fig. 6's micrograph shows it ending 5 rows short. The option
    `reproduction.pc_last_short` moves that end only in steps of about n − 1 rows, because each step adds a row to every
    middle leg. The Fig. 6 presets' 257 µm, the image's 5 rows, ends it 1 row short in the design preset (legs 1 to 7
    drawn 861 / 870 / 870 / 870 / 870 / 870 / 869 rows, so that its drawn end lies 4 rows beyond the image's) and 3 rows
    short in the as-fabricated one (legs 1 to 7 drawn 857 / 866 / 866 / 866 / 866 / 866 / 863 rows). The last leg's end
    changes no row's physics (each row is drawn where the schedule puts it on the path). It is a layout quantity that
    Gifford does not print, and the trace registers the last leg on its own along the flow (`provenance/rc6_traces/`).
26. There is no particle-passage check (rc6, raised on the Gifford Fig. 4 design-level zooms). opencif checks what a
    process can make (the floor, the corner limit, and the aspect ratio) and what the network extracts (the cutoff at
    every gap), and it does not check whether a particle that the gaps pass fits the side channel that it enters. A
    sealed start begins the side channel at the floor, which can be narrower than such a particle. For example,
    `gifford2014_fig4_design` draws 4.86 µm at row 1 against D 5.20 µm (its panels B and C show the pointed starts), and
    `gifford2014_fig6_design` draws 4.69 µm against D 3.62 µm (above it). Both are design-level pictures of the
    equations, already reported above the aspect-ratio limit (item 7), and their as-fabricated presets start at 20, 8,
    and 15 µm. Over the 46 shipped records that draw a device (`provenance/rc6_scripts/passage_scan.py`, which compares
    the narrowest straight side channel with D), 2 draw one narrower than D (`gifford_physics_linear`, 5.35 µm against D
    6.00 µm; `gifford2014_fig4_design`, 4.86 µm against D 5.20 µm), and 9 more draw one between D and 2 D
    (`m15_linear_shear_nomeander_filter_teardrop` 4.37 µm against D 4.00 µm;
    `m16_serpentine_shear_nomeander_plain_separate` 4.37 µm against D 4.00 µm;
    `m17_paperclip_shear_nomeander_filter_ring` 4.37 µm against D 4.00 µm; `m13_paperclip_const_nomeander_sealed_fine`
    4.37 µm against D 4.00 µm; `gifford2014_fig6_design` 4.69 µm against D 3.62 µm;
    `gifford2014_fig5D_as_fabricated_12p5um` 12.50 µm against D 9.03 µm; `gifford2014_fig4_as_fabricated_8um` 8.00 µm
    against D 5.20 µm; `xia2016_rs7_specified` 11.36 µm against D 6.86 µm; `gifford2014_fig5C_as_fabricated_12p5um`
    12.50 µm against D 6.98 µm; the scan reads each design's planned closure widths). The report of rc7 on the drawn
    mask's measured widths (`validation.passage`, described below) agrees on the two below D and, over the 50 rc7
    records, flags nine between D and 2 D. These nine are the list above without `xia2016_rs7_specified`, whose 11.36 µm
    is the closure width of row 3, a rectilinear transition row that the mask draws at the transition's own width (its
    narrowest drawn straight side channel is 16.02 µm, 2.33 D), and with `gifford2014_fig4_aligned` (8.00 µm against D
    5.20 µm, 1.54 D, and below the 10.23 µm that its rows 1 to 3 pass). The future work is a passage rule (the side
    channel at least a stated multiple of D wherever the rows extract), reported like the aspect ratio. Since rc7
    (ADJ-C9), every record carries the report (`validation.passage`, which compares each measured channel class with D
    and 2 D) and a warning note below 2 D (`checks.passage_report`, on by default). Where rows are realized above D
    (declared or not), every class is also compared with the largest particle that a gap passes, min(D_realised_max, G),
    with its own warning (code review R3). Because the rows that feed each class are not traced, this comparison bounds
    every class (for example, the Fig. 4 aligned start's 8 µm side channel against rows 1 to 3's 10.2 µm). It is still a
    report rather than a rule, because there is no particle-wall model.

27. Dinh's head as drawn (rc7). The quarter turns are square with Fig. 2B's roundings (DINH-1), whereas the units' far
    U-turns stay semicircles (centerline radius s/2). No micrograph shows them, and Fig. 2A's inset, a drawing, shows
    them square. Drawn square, they would add two corner terms and 0.43 s of centerline per U-turn, which the
    exact-length solve would absorb by moving each unit's reach by a few µm, and the resistance that each unit draws is
    unchanged. Every corner of the device on Fig. 2B panel 1 is rounded, as on the mask. The corners where the units'
    stems leave the trunk (the head gaps' entrance) measure about 7 to 15 µm against the 5 µm head fillet (11 µm drawn,
    `meander.junction_radius_um`; the reading moves by a few µm with the registration, and on the final rc7 trace the
    head's 12 wall corners read a median of 14.2 µm, from 4.4 to 31.5, the largest beside the panel's arrows and label,
    so the drawn 11 µm lies below the median). These corners are drawn only, because they enter the model through the
    gap entrance coefficient c_ent, an external input (1.0; the design is "insensitive while M13 passes", spec v5 S:756,
    and the record's "mouth R_gap ±20 %" scenario, item 9, is how it is represented), rather than through the T-junction
    term, whose corners lie at the wall's far face (register DINH-7). Those T-junction corners, which the term does
    model through its fillet part, read 7.6 µm on the panel (one clear reading, and another, near the panel's arrows, of
    18 to 33 µm) against the 5 µm fillet drawn and modeled. At 7.6 µm, the fillet part would move c_T by about −0.075 w
    (−2.1 µm of length per unit), about 1.3 times the term's band. In rc7, the corners are neither drawn nor modeled at
    7.6 µm, because that value is one noisy reading and the fillet sets every head corner, and the reading is disclosed.
    Two head readings are shipped, the printed f* (`dinh2024_v1_depicted`) and a two-cutoff reading of Fig. 2A's head
    and the printed length (`dinh2024_v1_figure_aligned`, register INV-1), whose head rows at 2.17 µm are declared
    departures below the κ plateau (the record notes this in `validation.gates.plateau_u_min_all_planned_rows`).
28. Gifford Fig. 4's layout (RC8; register GIF4-1, INV-2). The reproduction is the nested reading with the 8 µm start
    and the inset's measured gap (`gifford2014_fig4_aligned`, RC8). In this reading, the four legs are nested, as the
    caption ("the second and fourth legs"), the schematic, and the same paper's Fig. 6 say, as the inset's two arrows
    suggest, and as the inset's own taper direction confirms. The arrows point the same way, which is the weakest line
    of evidence, because their placement, symmetric about leg 4, also allows both to mark leg 4. The taper direction is
    the check V13. In a CIF leg, the side channels widen along the flow, so the sign of a leg's outer width slope is its
    flow direction. Through the trace's registration, the neighboring leg's slope on the image has the nested leg 2's
    sign (PASS; the image's slope is −0.00240 ± 0.00059 µm/µm, 4.1 SE, against the mask's −0.00257, at 130 of 130
    eligible positions with a ±30 µm window) and the opposite of reading B's leg 3 (FAIL; −0.00254 ± 0.00062, 4.1 SE,
    against +0.00301, at 117 of 117). Leg 4's own array part, the control, passes (PASS; −0.00469 ± 0.00040, 11.6 SE,
    against −0.00800, at 80 of 80 with ±12 µm). The check is a sign test whose standard error is the larger of a
    Newey-West one, with lag 100 µm, and the independent-sample one (`provenance/rc7_traces/code/taper_sign.py`). Every
    printed number is kept, and the frozen physics is exact on the drawing, with a declared off-design entrance. That
    entrance is Gifford's start rule at 8 µm, which is App. A's "other start mode" (App. A, lines 1118 and 1119) with a
    mass defect of 0.0022. Rows 1 to 3 lie above D, row 1 at 2.0 D, with 0.5% of F, and the 8 µm row-1 side channel is
    narrower than the 10.2 µm particles that those rows pass.

    One feature is not reproduced. Leg 2's side channels beside the inset read 17.5 to 18.1 µm per side wider on the
    image than the printed single-f schedule draws them. The schedule draws leg 2's rows there (1,730 to 1,758) at 105.1
    to 106.7 µm, and the image's 122.6 to 124.8 µm (the drawn widths plus the walls' measured offsets; the width
    profiles read 124.5 / 125.0 µm) are what it draws about 280 rows (12.7 to 13.0 mm of leg) further on, in leg 3
    (1.17 × the schedule's width). A registration or wall-finding error cannot be that large, because the pitch and the
    legs' spacing read to 0.1%, leg 4's walls read within one image pixel, and leg 2's posts and center line sit on the
    image's. No reading with the printed single f and a start inside the beads' bracket explains it, so the question is
    open. Reading B (the legs in path order, rc7's aligned preset) reproduces those widths only because leg 4's early
    end was fitted to them, and it is rejected by the taper direction, the caption, and the schematic (and by the
    arrows, read as one per leg). It is kept as the variant
    `provenance/rc7_traces/variants/gifford2014_fig4_pathorder.yaml`, a figure. The sealed design start
    (`gifford2014_fig4_design`, the design-level reference) would put every row at D, although its 4.86 µm first side
    channel at 80 µm depth is 16.5 : 1, against the "10 : 1 aspect ratio limitation" that Gifford prints (p. 4503). The
    8 µm is an inference (A), because Gifford prints neither his first side width nor a feature floor. The lab's later
    guidelines (Lezzar et al. 2022, which call features under 20 µm "very difficult to mass produce" and aspect ratios
    above about 7 : 1 "very limited") would give 20 µm at 80 µm depth by the feature size (11.4 µm by the aspect ratio,
    inside the beads' bracket below). However, the inset's own gaps read 17.2 µm, so 20 µm was not this device's floor.
    The inset's beads bracket the start (register INV-2; `provenance/rc7_scripts/postfreeze/early_loss.py`, a streamline
    model with steric exclusion). A straight start of up to about 12 µm is allowed (the chance of the 2 or fewer beads
    seen in the outer band is 0.97 at 8 µm and 0.095 at 12 µm), and starts of 17.2 µm and wider are rejected (2.8 × 10⁻⁶
    at 17.2 µm, 8 × 10⁻⁹ at 19 µm). A device held to a 20 µm floor can still start sealed, with its rows narrower than
    20 µm folded into a meander head (spec v5's FOLD start, Xia's remedy). The variant
    `provenance/rc7_traces/variants/gifford2014_fig4_fold20.yaml` does this, with legs 1 and 2 at the generator's own
    spacing for the head and legs 2 to 4 on reading B's rows (a path-order figure). Built at the RC8 code, it exited 0
    with N 3,546 (a 59-row head) and every check passing, no row's cutoff above the 8.3 µm beads, and, on the inset, a
    claim median of 1.10 µm and a p90 of 2.51 µm at its 20.5 µm gap, as the design preset's. It is shown as a figure
    rather than a reproduction, because Gifford drew no head.
29. Parts that collide (rc7 final; register BUILD-2). Before a photomask is written, the build draws the outline in
    memory and refuses it (exit 5, with the two parts, the shortfall, and the remedies, and, when figures are drawn,
    `<name>_refused_layout.png`) in any of the following cases:
    - A part (the inlet with the head, a leg, a U-turn, or the outlet) lies closer than the floor to a part that it is
      not joined to.
    - The outlet ports are not separated (from each other and from the inlet, in the drawn fluid).
    - An outlet port lies within 2 mm of the inlet port.

    `opencif check` and `generate(check_only=True)` draw the same outline and refuse the same way (this takes minutes
    with a head, and `check_only='plan'` plans without drawing). Where only the outlet meets the inlet end (a folded
    layout whose last leg ends beside its inlet, or a straight device too short for its ports), the build lengthens the
    outlet stretch by the ports' shortfall plus 200 µm, at least 500 µm a step, up to 12 times, and re-plans. The
    closure takes the longer stretch, and the record names it (`design.outlet_stretch_extended_um`). The build does so
    in every gate mode, reproduce included, except for a paperclip's ring outlet (the ring sets that stretch) or when
    another pair of parts collides too. It stops when a longer stretch does not help or no longer closes the outlet, and
    the gate then refuses, with the refusal listing the tries (`refused.layout.outlet_stretch_extensions`). A layout
    value that the spec gives (a leg spacing, measured leg axes, or an inlet offset) is used as given and refused if it
    collides, and it is never moved. With `legs: auto`, the fitter skips the counts whose drawn parts collide, and when
    every count collides, it hands the build the first count whose outlet alone meets the inlet end. It judges counts
    without the build's floor and fail-safe re-plans, so a count that it accepts can still be refused. What the
    collision check does not cover is as follows.
    - Parts joined in flow order (the head and leg 1, or a leg and its U-turn) are not distance-checked against each
      other, and nothing checks a contact between them away from their joint as such. Such a contact may show on the
      drawn mask as a changed hole count, a lost wall, or a changed measurement, or it may not show at all.
    - Those changes, and every other integrity check of the drawn mask (side-channel continuity, every planned feature
      found, the corner radius, and the DRC), FAIL the build, and a FAIL's photomask is written as `<name>_FAILED.dxf`
      and never under the deliverable name (the build writes `<name>_UNASSESSED.dxf` until its verdict exists).
    - For the ports, the 2 mm rule applies between the inlet port and each outlet port. The outlet ports are checked to
      be separated from each other rather than for a distance, and the inlet port is kept the floor, rather than 2 mm,
      from the other parts.
    - For arrays, before the array's photomask is written, the copies are checked for merging or disconnection (one
      component each, or one in all with a shared inlet), each copy's product port is checked for separation, and the
      manifold is checked for touching the copies only at their ports. RC8 (adjudication of the rc7 reviews, D4; K8,
      K11) applies the 2 mm rule between copies, from each copy's own inlet port as drawn (Dinh's filter inlet port has
      a 1,020 µm radius) or the shared inlet port to every other copy's outlet ports, gated with the other checks before
      the array's photomask is written (`checks.inter_copy_port_spacing` in the array's record,
      `compose.inter_copy_port_spacing`). The outlet ports of neighboring copies are reported there rather than gated,
      because no rule sets a distance between outlet ports within one element either (742 µm on both shipped elements).
      The shipped arrays place them 730 µm (`m07…_x3`) and 300 µm (`m13…_x4`) apart, edge to edge, so tubing for the
      outlet ports of an array is the user's to place (rc7 had no rule between copies at all; code review 4, F10).
    - For contradictory inputs, `spec.conflicts_with_reproduction` refuses the spec keys mapped in `spec.SPEC_ENGINE`
      and the code review's pairs. RC8 adds two pairs (Review 2 of rc7; K7). The first is `meander.enabled` against a
      contrary `reproduction.head`, either way round (`enabled: false` beside a head that the reproduction draws, or
      `enabled: true` beside `head: none`). The second is an explicit width that sets the head's arms
      (`meander.narrow_um` for the narrow style, `meander.channel_um` for pack and arms) beside a different
      `reproduction.w_arm`, with one message per conflicting head. The second pair is a policy, one source for the
      head's widths, rather than a contradiction in the engine's terms, because `narrow_um` sets the N rows' width and
      `w_arm` sets the A rows' width and the arm pitch (`opencif/physics/schedule.py`, `Spec.w_arm` and `Spec.w_narrow`). Since
      `w_arm` defaults to the narrow width, the refusal stops an explicit width from silently differing from the arms
      that it would otherwise set (review 3, D3). A two-width head stays expressible in the reproduction block alone.
    - `tests/test_composition.py` samples the composition space (layouts, heads, inlets, outlets, start rules, four
      processes, and given and derived spacings) and builds a subset, so a combination outside the sample is covered by
      the refusal rather than shown to compose.

    Since RC8 (K10), each of the extension's tries re-plans from the spec's own `cap_frac` and side inlet width, and the
    kept try's values are restored with its stretch, so a build's photomask no longer depends on the tries before the
    kept one (rc7 re-planned from the values that the previous try left; Review 1 of rc7). A later try that can no
    longer meet the floor ends the loop, as an outlet that no longer closes does, with the refusal listing the tries.
30. What `check` predicts (RC8; K13, K14). `opencif check` and `generate(check_only=True)` answer whether a design is
    inside the validity envelope and whether its parts clear one another, and "inside the envelope" is not a prediction
    that the drawn mask passes Q2. Two plan-stage cases are known where the drawn mask does not pass Q2, and RC8 says so
    in the check's and the build's notes before anything is drawn:
    - The first is a floor start outside reproduce mode (`start: floor` in strict or evaluate mode, without a meander
      head). Its entrance rows are declared departures in reproduce mode only (`build.active_options`), so where they
      over-extract, as they usually do, the build FAILs Q2 on undeclared rows (the shipped `m14…_floorstart19`, in
      evaluate mode with exit 2, is such a FAIL by design). A floor start just below the floor may still pass. The
      sealed-start refusal's remedy names this case too.
    - The second is a plan whose filtrate fraction overshoots the target by more than `checks.F_tolerance`. The
      schedule's march stops at the first row whose cumulative fraction reaches F, so it can overshoot by up to one
      row's extraction (`schedule` "opencif F10"), and Q2's F check is then expected to FAIL. `physics.trim_F` trims the
      last gap to F exactly. No shipped plan overshoots by more than 0.00040 (`dc_D5_F0p3`), which is inside the default
      10⁻³.

    A design refused by the plateau gate in strict mode is refused before its head's sealing search
    (`pipeline.plateau_preflight`, RC8), because the gate needs only u = D / w_c0, which every drawn record equals.
31. Found by the independent review of the 1.0.0 code (2026-10-23) and left for the next release, because each change
    would change the package source and so the code stamp, which this release keeps. None of them moves a shipped
    photomask, record, or verdict.
    - The width inventory of the aspect-ratio check, P9, and P9-AR takes the head units at their planned widths and
      leaves out the head's gap mouths (item 7); the fix takes the measured unit widths (`measure`'s per-unit minimum)
      and the head-mouth gaps.
    - Strict mode refuses a design whose plan is outside the envelope, but an envelope check that fails only on the
      drawn mask (the physics audit's M14, or a gate recomputed on the drawn plan) gives exit 3 in strict mode too. On
      the shipped strict designs the drawn lamina's Reynolds number is at most 1.0031 times the plan's, so only a plan
      within about 0.3% of the Stokes gate is exposed. No test covers this path.
    - A head unit's measurement (item 12) should integrate r(w) along the unit's centerline, as the U-turns are measured.
    - The DRC passes compact solid islands narrower than half the floor (4 to 9 µm squares at a 19 µm floor are not
      flagged, 10 µm squares are); the topology counts catch an extra island, and R6b checks the straight rows' posts,
      but the sizes of the filter's and the outlet's islands are not checked.
    - The constant-shear taper's table spans [max(`width_min_um`, 0.3 × `width_in_um`), `width_in_um`] and clamps to
      its ends without a note, so a `width_min_um` below 0.3 × `width_in_um` gives a floor of 0.3 × `width_in_um`, and
      a `width_min_um` above the inlet width is accepted (strict mode's M16 refuses all but a sliver of that case).
    - Untested paths: strict mode's exit 3 above; a build with the DRC off at the pipeline level (exit 4); Q1 on a
      drawing defect seeded inside a build; and the standalone audit with a DRC that errors, which gives FAIL (exit 2)
      where a build gives INCOMPLETE (exit 4). The fast suite that pins the closure checks only constant-width designs
      without a head; the closure on a taper, a head, and U-turns holds (to 4 × 10⁻¹¹ on m15 and m17) but is pinned only
      by the long suites' Q1 and Q2.
    - The drawn-mask M14 compares values rounded to four significant figures, so it passes up to 5 × 10⁻⁵ relative
      above its limit, whereas the plan's gate compares exact values.
    - The outlet leads' two end zones (1.5 w each) are counted at their design width, and the stretch's side width is
      one section at its midpoint, although `measure`'s docstring says the leads are measured with no trust in the
      design (item 15a gives the last rows' sensitivity).
    - `opencif build` writes its summary in the system's default text encoding, which on Windows is not UTF-8, so a
      build there writes its photomask and record and then stops with exit 1 at the summary (it contains δ and κ). On
      Windows, run Python 1.0.0 in its UTF-8 mode: `python -X utf8 -m opencif ...` (README, "Installation"). The
      photomask is written in text mode, so on Windows its lines end in CR LF and its sha256 differs from the
      release's, while its geometry is the same.
    - `schedule.turn_side_width` is unused and would take the largest U-turn radius, and `bend_ratio` raises a radius
      ratio below 0.01 silently (reached only by the root finder's probes, never by a drawn arc).
    - Outside the package, `tests/test_trace_anchors.py`'s docstring still says that the publishers' figures are not
      in the repository. Since 1.0.0 they are in `notebooks/figures/` (`NOTICE.md`); the docstring is left as it is so
      that the tests keep their stamp.

## E. Known cosmetic issues, deferred past the freeze

Each issue would change the package source and so the code hash, and it would need a full rebuild of every record. None
of them affects a design, a check, or a verdict. (Release rc2 fixed three earlier items, which were the absolute paths
in DRC audits, the plateau band's `method` text, and the nested parentheses in the band text.)
- The default name `cif_D<D>_F<F>` repeats for two designs with the same D and F, so pass `name=` (or `--name`) for such
  designs.
- Separate-outlet records carry an unused `design.spec.outlet_style: 'teardrop'`, which is kept because the replica
  factories that `tests/test_presets_match_replicas.py` compares against carry the same default.
- The command line prints a NumPy RuntimeWarning (from `opencif.physics.kernels.r_inverse`) before "SPEC ERROR" when a closure has no
  solution.
- Some zoom markers on the layout figure are drawn as a single edge when the zoom window is taller than the strip.
- (rc4) The pre-drawing w_min check reads the planned side widths, so with the papers' rectilinear transition
  (`meander.transition: paper`, a declared departure) it can name a planned width that the drawing replaces. For
  example, `xia2016_rs7_specified` notes "straight side channel w_s = 11.32 um (row 3) below w_min ... 16.00 um",
  whereas its drawn narrowest straight side channel is 16.0 µm and the drawn aspect-ratio check passes.
  `transition: paper` is a reproduction-only option (in `spec.py`, it needs `gates.mode: reproduce`), so the note can
  appear only in reproduce mode, where it is a note and never a refusal.
- (rc6) A build whose outlet stretch ends beyond the solved 2-D tables says so in its verdict only when it passes or
  reproduces (exit 0 or 3, as with the other 1-D model notes), whereas a FAIL or INCOMPLETE verdict names its own
  causes. Every record carries the status (`design.outlet.stretch_end.table_domain`). No shipped design reaches the
  tables' ends (the largest width ratio is 3.87, on `wbc_paperclip_teardrop`, against the steps' 8 and the flares' 4).
- (rc6, from the second independent review; fixed in rc7, ADJ-C8) A spec that sets a value through a reproduction block
  (for example, `reproduction.r_fillet`) listed the spec key's default or derivation among the record's `defaults_used`.
  Release rc7 names the reproduction value as the one drawn and gives each derived value its effective value.

(RC8 fixed the two rc7-final items. `pipeline.name_mask` now names the file just assessed and removes the other
verdict's stale file, and `compose.array` clears the array's own outputs before composing, so an earlier photomask can
no longer stand under the deliverable name or replace a new FAIL. Every full build records `validation.drc.file_final`,
the name that the verdict gave the photomask, beside `validation.drc.file`, the name that the DRC audited
(`<name>_UNASSESSED.dxf`; the same bytes, `files.dxf_sha256`). RC8 also refuses, with the field named, input that rc7
let through (rc7 refused a negative or zero depth, although only from inside the planner and without naming the field).
The refused inputs are a negative, zero, infinite, or true/false length, depth, flow, radius, tolerance, or count (37
keys must be > 0, and the 9 that may be 0 are the fillet, the gap-entrance coefficient, the chip margin, the three
corner and junction radii, the inlet port gap, the flare start, and the tubing length), including `chip: (0, 0)` and
`array.bus_width_um: 0`, which rc7 read as "no chip" and "the default bus", `generate(array=True)` (one copy in rc7),
and a non-finite D. An unknown process and `seal` without a head, which raised `ValueError`, are `SpecError`s (review 3
of the RC8 candidate, E9; Review 1 of rc7; patch review n1).)

(The final freeze of rc7 fixed an item that `tests/test_composition.py` found. Where a leg's last wall block runs past
the leg's end into a U-turn (a short paperclip whose leg 1 holds only the head, for example, the sweep's Gifford-process
path-order case at F 0.5), the arcs of its upstream corners were cut at the leg's end, leaving a 0.19 µm step against
the turn's wall that the corner-radius check FAILed. The arcs now end on the band edges before the leg's end
(`layout.leg_linear_local`, register BUILD-2). No shipped device ends a leg on a wall-block row, so none changes.)

(Release rc5 fixed the rc4 item on two stale code comments, the `api.py` module docstring and `pipeline.wmin_check`'s,
which described the drawn aspect ratio as reported only. It also writes G9's leg path into the record as
`validation.length_budget.definition`, with leg 1 running from the head's top and the last leg to the end of the outlet
stretch, so the reported path is at least spec v5's ΣΛ, as a review of rc4 noted.)
