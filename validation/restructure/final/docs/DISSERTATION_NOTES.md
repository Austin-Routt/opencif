# What the dissertation must say about opencif

This note maps the code base onto the dissertation's frozen foundation, which consists of Appendix A (step-by-step
derivations), Appendix B (the argument table, claims C1 to C9, and the thesis), the glossary
(`reference_glossary/PREAMBLE_FINAL.tex` holds its numeric macros), and the design specification
`Ladder_Array_Design_Spec_v5.md` ("spec v5"). All four, with their md5, are in `foundation/`. For each idea, the note
gives the source to cite and the code that implements it, and it then states what is new, what may be claimed, and the
notation map. The statement-by-statement traceability is in `docs/ALIGNMENT.md`, the limitations and future work are in
`docs/LIMITATIONS.md`, the results are in `docs/RESULTS.md`, the review history is in `docs/REVIEW_2026-10-08.md`, and
the source of every preset value is in `docs/PRESET_PROVENANCE.md`.

The note follows four citation conventions. App. A's equations are numbered as printed in `APPENDIX_A_revised.pdf` (R =
Result). App. B is cited by claim and premise (e.g. P-C7.1) or by critical question (CQ), and App. B's own "B5"/"B6" are
regression stamps of its reference script rather than sections, so they are never used to cite opencif. Glossary
subsections are cited by name (e.g. "glossary, Validity-Envelope and External Inputs"), because "A.n" inside the
glossary prints as "C.n" in the assembled dissertation. Finally, opencif is the dissertation's Prospectus C7 (the
targets-to-photomask tool) and C8 (the machine-enforced envelope), so never call it App. B node C7 or C8.

## 1. Where the code follows the frozen framework

| idea | cite | code |
|---|---|---|
| the design targets are the size cutoff D and the filtrate fraction F (in App. B, "yield" means F, the fraction of fluid sent to the filtrate, and never particle recovery); the other inputs are the channel shape, the depth, and the gap and post sizes, whereas f_gap and N are outputs | App. B thesis; App. A R13 (setup) and R13(c), eq. (86) | `spec.py` (`targets.cutoff_um` and `filtrate_fraction`; optionally, the flow ratio F/(1−F), which is not the glossary's FR, Inglis and Herman's filtration ratio); `api.generate(D, F, ...)`. Only reproduce mode may prescribe a paper's depicted head fractions (`reproduction.f_head`, `depicted_head`), and those rows are declared |
| one declared kernel pair, used at all three call sites, which is the cubic f = 3x² − 2x³ with κ_g = 4.9138 [4.750, 5.078] (half-width 3.4%) and never κ_F | App. A R1 (eq. 10), R2 (eqs. 18 and 20), and the R9 corollary (eq. 59); glossary, Cutoff-Law Prefactor Decomposition, pairing rule; App. B C8 kernel declaration | `opencif/physics/kernels.py` (the design code's kernels, restated and tested bit for bit against the frozen `reference/CIF_v2.py`), `schedule.f_of`, `x_of`, and `KAPPA_G`; every record carries `design.kernel_pair` |
| the lamina target δ\* = D/κ_g; every gap extracts f(i) = g(δ\*/w_c(i)), so the lamina is δ\* at every row, and on the κ plateau every designed gap's cutoff is D | App. A R13(a) (eq. 81), implementers step 2; App. B P-C8.1 | `spec.to_engine`, `schedule.build`; `tests/test_appendix_trace.py` |
| D read on the exact cubic root (eq. 20), as App. A R13(a) does, whereas κ_g was fitted on the leading-order δ = w_c√(f/3) (eq. 18); the two differ by at most 0.47% up to κ_g's calibration cap f = 5.76 × 10⁻⁴ | App. A R2; glossary, Equations: Step-by-Step Unit Derivations, "Skim cubic and dividing-streamline width" | `schedule.x_of` (tested) |
| the march P(i) = (1 − 2f)P(i−1), S(i) = S(i−1) + f P(i−1), and the mass identity P + 2S = 1 + 2σ₀ | App. A R10 (eqs. 60 to 67) | `schedule.build`, tested (asserted < 10⁻¹², observed 4 × 10⁻¹⁴). A sealed start stops on F = 1 − P(N) (eq. 67), whereas a floor start stops on the F that its drawing delivers, 1 − P(N)/(1 + 2σ₀) (§2.4) |
| D and F are independent inputs on the plateau, because D fixes δ\* and F then fixes N (a lower-triangular Jacobian); N·D² is invariant at fixed F | App. A R13(b) and (c) (eqs. 82 to 86); App. B C8 (P-C8.2) | `schedule.plan`. The tests show that changing D changes N and not F, that changing F changes N and not the cutoff, and that N·D² is constant to 0.51% over D from 3 to 7 µm on w_c 100 µm (D 3 µm is just below the plateau, u 0.030) and lies within 0.92% of eq. (82). Do not quote App. B's 0.29% for opencif, because that is the F_lam pair's number |
| the transport lemma 2φ_crit = F | App. A R12 (eq. 80); App. B C9 | tested with `CIF_v2.critical_phi` (observed ≤ 1.4 × 10⁻¹⁴) |
| the rectangular-duct resistance (Gifford eqn 5) | App. A R14 (eq. 88) | `schedule.rr` → `kernels.r` (= `CIF_v2.r`, bit for bit; 50 odd terms, sorted), tested against a 401-term series and the square-duct constant 28.454154 |
| the validity envelope, made of the κ plateau u ≥ 0.035 and its foot 0.012 (external inputs; u is the glossary's reduced size ratio d̃ = D/w_c), channel Re ≤ 10, and the Stokes gates Re_gap ≤ 10⁻² and Re_δ ≤ 10⁻² | glossary, Validity-Envelope and External Inputs; App. B C8 CQ7; App. B thesis, "The Stokes regime is gated in code" (CIF_v2 M14; spec v5 §19 M14) | `build.gates`; strict mode refuses before drawing (`pipeline.prepare`, exit 5), and the refusal tests are in `tests/test_negative.py` |
| gap dominance M13, max R_gap f/R_c < 0.05 (not an App. B gate; the v114 glossary identifies it with the "perturbation" to gate on) | spec v5 §19 M13; CIF_v2 M13; glossary, entrance length | `build.gates` |
| the thin lamina δ(i)/h ≤ 0.10 at every planned row (E11, spec v5 "∀i"; rc4; CIF_v2 checks the cutoff rows only) and the taper rate \|w_c(i) − w_c(i−1)\|/(G + g_p) ≤ 10⁻² (M16) (not App. B gates; rc3) | spec v5 E11 and M16; CIF_v2 `march()` with `Tolerances.eps_lam` and `eps_taper` | `build.gates` `lamina_ok` (cutoff rows), `lamina_all_rows_ok`, and `taper_ok` (strict refuses); the trace shows parity with `CIF_v2.march` just below and above each limit |
| CIF_v2's floor w_min = max(F_min, h/AR_max) on the gap, a drawn GS, every center width, and the terminal side width (rc3), and spec v5 G7 on every straight side channel and the drawn mask (rc4; G11 would fold a narrower side channel, whereas in strict mode opencif refuses such a plan and fails (exit 2) a drawn mask over the limit) (not an App. B gate) | spec v5 §18, G7, and G11; CIF_v2 `Process.w_min`, `design()`, and validator G7; glossary (AR_max external) | `pipeline.wmin_check` (strict refuses, and other modes report); the `assess` aspect-ratio check (in strict mode a pass condition, and in other modes of kind `fabrication`, reported); LIMITATIONS 7 |
| the fold width over the gap, GS/G, in γ_GS [1.4, 1.5] when a meander head is drawn (not an App. B gate), and f_soft = 5 × 10⁻⁴, reported | spec v5 §2.3; §14.4 and §18 (pre-flight "ASSERT GS/G in gamma_GS"); §18 "if max(f) > f_soft: WARN" and §19 validator D3; CIF_v2 `design()` and `Pairing`; glossary, Validity-Envelope and External Inputs (both external) | `pipeline.gamma_gs_check` (strict refuses, and other modes report); `build.gates` `f_soft` and `f_above_soft` |
| the length budget is reported rather than enforced, so a design over its chip area is drawn, and the verdict says so (opencif reports the 2-D footprint against the chip and, since rc4, also spec v5 G9's own predicate, each leg's path length against the chip's leg length, `status.length_budget_G9`) | App. B thesis ("it reports, without refusing, a design over its length budget (C8, CQ5)"); C8 CQ5 ("reports the result without failing") | `pipeline.assess` (check kind `budget`; `assessment.over_chip_area`, `api.Result.over_chip_area`); `compose.array` |
| below the plateau, the cutoff is reported as App. B C7's band, with κ from the frozen Hill fit (κ_max 5.568, u_½ 0.01651, n 1.791) at u = κ_g δ\*/w_c and the band [κ(u) δ\*, κ_g δ\*]; for Dinh V1 this is [1.80, 2.17] µm, App. B's band for a lamina held at δ\* (its Xia-scaled reading), whereas App. B records the wider [1.25, 2.17] because Dinh's schedule is not printed. The band is plateau-anchored rather than the self-consistent cutoff, which can lie lower (App. B C7 CQ3), and App. B warrants it only as a low-severity consistency band | App. B P-C7.1, CQ3, and CQ4; glossary macros `\hillKmax`, `\hillUhalf`, `\hillN`, `\dinhCollapseLo`, and `\dinhBandLo/Hi` | `build.gates` (`cutoff_band_um`); `tests/test_appendix_trace.py` reproduces 1.25, 2.17, 4.08, and 1.80 |
| the generator stops with an error if δ\*/w_c ≥ ½ | App. A, "For implementers", step 3; spec v5 §19 D4 | `spec.to_engine` (SpecError), `schedule.f_of` |
| supplied constants that no validator can check, which are c_ent (= 1), R_min, σ₀, and (rc7) L_eq, here the square corner's computed c_corner | glossary, Validity-Envelope and External Inputs; spec v5 §20 | coded A (c_ent, R_min, σ₀) or G (c_corner, a 2-D computation) in `PRESET_PROVENANCE.md` |
| the segment rule, by which a gap-free span (a U-turn) carries its upstream rung's flows and takes that rung's ρ | spec v5 §1 and §14.5 | `schedule.reclose`, `layout.loop_profile` |

## 2. Where the code generalises the framework (cite both)

1. For the closure, App. A R10's closure (eq. 69, which with eq. 68 gives Gifford's eqn 4, eq. 64) assumes pressure
   equilibration, whereas opencif closes on the exact flow-prescribed pressure field with gap drops, the 'exact gap-drop
   closure'. It sets P_s(i) = P_c(i) − R_gap(i) Q_gap(i) and R_s(i) = (P_s(i) − P_s(i+1))/S(i), which is the identity
   R_c Q_c − R_s Q_s = R_gap Q_gap(i) − R_gap Q_gap(i+1) of spec v5 §13 (its equivalence chain) applied row by row. The
   closure reduces to eq. (69) when R_gap Q_gap is equal at consecutive rows, and it stays exact where the gap geometry
   or the row pitch changes (the head arms, the seam, and the U-turns). In a test, the network solved independently on
   the closure resistances gives every designed f (observed ≤ 5.7 × 10⁻¹⁰).
2. The outlet closure is exact as well. opencif does not use the outlet identity F = 2R_c/(R_s + 2R_c) (App. A R8 eqs.
   44 to 47; R11 eq. 76, which needs the PE outlet closure). It sizes the outlet leads from the last row's exact
   pressures instead, so the network delivers the design F (observed ≤ 6.4 × 10⁻¹¹). The march itself satisfies eq.
   (75), 1 − P(N) = 2σ(N)/(1 + 2σ(N)), to rounding, whereas eq. (44)/(76) evaluated on opencif's own outlet resistances
   differs from F by the outlet row's gap drop, about 3 f_gap, so the implementers' step-7 outlet check does not apply.
   App. A R8's last paragraph (P3 holds only approximately in a finite network) and App. B C4 P-C4.3 (network F equals
   the identity to 0.12%) are the context. Do not call opencif's outlet "the outlet identity".
3. opencif evaluates the Stokes gates at every row. CIF_v2 M14 evaluates Re_gap and Re_δ at row 1, whereas opencif
   evaluates both at every row at the design values (f\*(i), δ\*) and includes CIF_v2's row-1 term, so its gate is never
   weaker. Re_gap is the gap flow's Reynolds number in the gap section (CIF_v2's `reynolds`). Its ratio to the flux
   product f·Re_channel that App. B's thesis section and the glossary quote is (w_c + h)/(G + h), 1.2 to 2.8 on the
   presets, which makes it the stricter reading. Rows planned away from f\* (an unsealed gap 1) claim no cutoff and are
   reported rather than gated, as CIF_v2 does for prescribed schedules.
4. For the floor start, App. A eq. (67) gives the schedule's F = 1 − P(N) for any σ₀. opencif's floor start draws a side
   channel that starts at the floor width but carries no flow at its sealed end, so the drawing delivers 1 − P(N)/(1 +
   2σ₀). The march stops on that value (tested; Gifford Fig. 4 at 20 µm then ends at 299.9 µm against the printed
   300 µm), and the entrance rows that the floor over-extracts are predicted and declared. App. A's eq. (67) is the
   schedule's bookkeeping and holds. Spec v5's revision 5 called 1 − P(N) the physical F for any σ₀, which a sealed
   drawing does not deliver; its revision 5.1 (2026-10-23) states the delivered fraction and stops on it, as opencif
   does (erratum D108, `docs/ALIGNMENT.md` S3). The design code's docstring still says "for ANY sigma_0".
5. The retention numbers of Gifford's Fig. 5 devices quoted below are opencif's own reading of published data, because
   they come from this project's reading of Gifford's printed enrichments (`provenance/gifford_fig5_d50.py`). They are
   new evidence and are not in App. A, App. B, or the glossary (whose Selectivity-Axis Validation uses Dinh and Xia and
   finds neither a sharp held-out test).

## 3. What is new and must be written up (Ch. 6, with the evidence named)

1. The curved-duct resistance is new work, which is not in App. A or B. It is the curved-coordinate counterpart of App.
   A R14, solved by the same depth-wise sine expansion as R9/R14. A U-turn, a meander hairpin, or an outlet arc is not a
   straight duct of the same centerline length. At Re → 0, the flow in a rectangular duct curved in its plane is purely
   azimuthal and solves ∂²u/∂r² + (1/r)∂u/∂r − u/r² + ∂²u/∂z² = −1/r, and a sin(nπz/h) expansion gives a modified-Bessel
   problem per mode and the flux in closed form (`schedule.bend_ratio`). We tested it against the planar-annulus (h →
   ∞), Hele-Shaw (h → 0), and straight limits (`tests/test_units.py`) and against an independent second-order
   finite-difference solve of the same equation (`tests/test_bend_fd.py`, within 3 × 10⁻⁷), and reviewer R3's own solve
   agreed to 2.5 × 10⁻⁸. The size of the term shows when it is left out, because without it Gifford Fig. 4's rows after
   the first U-turn move by up to 5.2% of D, and the most sensitive rows of the headed designs (at the seam, or next to
   a U-turn) move by 13 to 78% of D.
2. The junction terms are fits and are opencif's own. A bend's two straight-arc junctions remove κ_j·w of the arc from
   the developed correction (κ_j 0.50 deep … 0.33 shallow), and a head T-junction is shorter than the centerline
   convention by 0.28 to 0.47 w for sharp corners, plus a fillet term. Both are fits to 2-D Stokes solves in two limits,
   placed between them by a depth weight. They are the weakest physics in the pipeline, and they matter most at a head's
   seam rows (without the T term, m09's seam rows move by up to 65% of D). Neither appendix discusses them. App. B C4
   CQ1's "sibling" is a different junction effect, the 3-D gap-mouth effect on the per-gap f_eff, which is also
   unmodeled here (LIMITATIONS 9). Whether such residuals share a common 3-D origin is an open question there and in the
   glossary rather than a claim. Every record reports the rows that these terms' range can move (the model band). The
   scripts and logs of the 2-D fits are in `provenance/junction_fits/`. The fillet (`process.fillet_um`; rc6 default
   min(5, 0.4 × the narrowest drawn feature, 0.3 × the narrowest head unit), 5 µm on every shipped spec that uses the
   default) is drawn over the inlet-and-head zone and the outlet zone only (the legs' posts are stadiums, and the
   U-turns are arcs). A fillet at or above half the narrowest feature is refused, since the smoothing erases channels no
   wider than twice the fillet. `fillet_um: 0` draws the papers' square CAD corners and builds. In a review sweep on
   rc3, a serpentine D 4 µm / F 0.5 device on the generic 100 µm process and a 25 × 75 mm chip has c_T −0.437 w and a
   seam band of 16.4% of D, whereas the same device with square corners has the sharp value −0.293 w and 9.6%. Dinh's
   fabricated head (Fig. 2B, first panel) shows its corners rounded, although his CAD inset (Fig. 2A) is square. rc7
   reads the quarter turns' outer corners at 6.0 / 8.6 µm and their inner (wall-end) corners at 9 to 11 µm, the head
   gaps' entrance corners at about 7 to 15 µm (on the final rc7 trace, the head's 12 wall corners read a median of
   14.2 µm against the 11 µm drawn, so the drawn value lies below the median), and one T-junction corner at 7.6 µm
   (registers DINH-1 and DINH-7; LIMITATIONS 8 and 27).
3. The drawn mask is validated by two questions. In Q1, the 1-D network solved on resistances measured on the drawn
   polygon matches the intended network (each row within 0.25% of D in cutoff). In Q2, every gap of the drawn device
   realizes D ± 1% and F within 10⁻³. A row is exempt only when a declared option is shown, by intended networks with
   and without it, to move it. Beside them run checks whose expectations come from the request rather than from the
   generator. These require no intersections (every part the floor from every unrelated part, and the outlet ports 2 mm
   from the inlet port) and the drawn components as requested (head style, sealing, filter posts counted on the drawing,
   outlet, and legs). Since rc7's final freeze, the clearance checks run on the outline in memory before the photomask
   is written, so parts that collide are refused and never written, and `opencif check` refuses the same way (register
   BUILD-2). Claim it as a gate on every build, with `tests/test_composition.py` as a sample of the composition space
   (22 devices on four processes, 96 checks) rather than as a proof over every combination (LIMITATIONS 29). The
   evidence that the checks can fail is `tests/test_negative.py`, which seeds defects (with their evaluate-mode and
   gate-function follow-ups; rc5 adds the flare's spec errors) and synthetic clearance violations, and every one is
   caught (63 checks with the controls in rc7). The strict refusals among them are the κ plateau, channel Re, the lamina
   Stokes gate, γ_GS, E11, M16 (with M13 and the Stokes gates as a design, and alone on the gate function), w_min on the
   gap, and, since rc4, G7 on a straight side channel. In every mode, they also include a head that cannot seal gap 1
   under `seal: require` and, since rc4, a ring outlet that cannot close the design. `tests/test_units.py` shows that
   the Q1/Q2 metrics catch a 5% f error and a reverse-flow row on synthetic networks, and among the records, m14 is a
   real Q2 FAIL and wbc_xia_process_refused a real M13 refusal. The DRC, topology, side-channel continuity, corner
   radius, and head envelope checks, and the Re_gap-only and foot refusals, are not seeded, although their failures were
   seen on real drawings during development (review rows 38, 40 to 43, and 46). The chip area is reported and is never a
   failure (a control in the same suite).
4. A single generator call builds a device. `opencif.generate(D, F, **parts)` (and `python -m opencif generate D F ...`)
   builds the spec, runs every gate and validator, and returns the verdict, the DXF, and the figures, and
   `docs/OPTIONS.md` matches every part style to its code (generated from the schema and tested to be complete).
5. Also new are the resistance-exact wall rule (the straight side wall drawn so that its integrated resistance is
   exactly the closure's, whereas the earlier linear rule left errors up to 15% in f), the packing arms head (the
   default, which is Dinh's topology with as many arms as the width holds, at the narrowest width that seals gap 1), and
   parallel arrays of a validated element with per-element outlets.
6. In rc6, we traced the generator's masks on the papers' micrographs (`provenance/rc6_traces/`,
   `tests/test_trace_anchors.py`). Each panel that shows post rows, head units, or loops is traced with the shipped
   mask. The scale comes from the panel's own bar read with its anti-aliased ends (or, where no bar sets it, from the
   panel's printed widths), and the rotation is measured on the image. The translation is fitted twice, on features
   indexed from a non-periodic feature and independently by chamfer matching. Along a periodic array, the chamfer scores
   one pitch away about as well, so the along position is the non-periodic feature's, which Fig. 6's END_CHECK verifies
   row by row. Each panel is then judged against nulls (the mask shifted half a pitch along and 5 µm across) and a scale
   scan. The traces show Dinh Fig. 2B's four panels at medians of 0.44 to 1.29 µm, with the outlet's widths within
   0.5 µm of the print and the pitch the image's. They show Gifford's Fig. 4 inset at 1.21 µm, and Gifford's Fig. 6 legs
   identified from their side widths and registered at 1.4 to 1.9 µm each and, with their turns, on one map at 1.47 µm.
   They also show Xia's Fig. 3 panels c to g at 2.3 to 3.8 µm (relief contrast), and Xia's whole device on Fig. 3a,
   where the shipped schedule is 5.4 mm short and the variant with the cutoff that the drawing implies (D 6.5315 µm)
   fills it. Write this up as a test of the generator given its inputs, with the inputs' provenance stated. The inputs
   are the papers' printed values; values measured on the same figures (the pitch and post sizes, the layouts, the
   widths that set a panel's scale or place it along the device, and Dinh's flare), whose agreement is consistency
   rather than a test; and the generator's own choices, drawn in the no-claim color. The limits of `trace_anchors.json`
   were set with the results in view (regression limits), whereas each panel's scale and rotation and the agreement of
   the two registrations across the flow are blind. The presets follow one convention (PRESET_PROVENANCE), in which
   printed values are taken as printed and unprinted ones are measured on the figures where the figures show them. Thus,
   the Gifford presets take their figures' layouts (leg axes and spacing) and Xia's take Fig. 3a's, while Xia's printed
   targets stay and its figure-read targets are traced variants (D 6.5315 µm; F 0.94072 with D 6.5705 µm). Inputs read
   on a figure reproduce it and are not independent tests of themselves. The medians above are rc6's, on rc6's masks
   and window. Quote the shipped masks' traces instead (rc7 and RC8, the window rule of `provenance/rc7_traces/`,
   whose README tabulates them): Dinh Fig. 2B at 0.44 to 1.40 µm on both Dinh masks (the head panel's metric null
   ratio is 1.17, below the 1.3 criterion), Gifford's Fig. 4 inset (leg 4) at 0.61 µm on the nested
   `gifford2014_fig4_aligned` and 1.10 µm on the design preset, Gifford's Fig. 6 legs at 0.51 to 0.64 µm on the aligned
   preset and 1.09 to 1.87 µm on the design preset, and Xia's claimed panels at 2.02 to 2.78 µm (the others are
   recorded, not claimed). The candidate central figure is
   `figures/fig_dinh2024_fig2B_traced_rc6.pdf`. Gifford's Fig. 6 (`figures/panels/gifford_fig6_four_legs__rc6`, the four
   legs and their turns on one map, and `figures/fig_rc6_fix_gifford6_turns`) and Xia's Fig. 3a whole-device
   registration (shipped and variants) show the layouts. They are read on those same figures, so they are consistency
   rather than tests, except where along the flow Fig. 6's turns start (5.7 µm from the image's). `full_device/` shows
   each generated device whole, with the micrographs placed by their traces. Dinh's and Gifford's figures need RSC
   permission (README there), whereas Xia is CC BY 4.0. The author's review of 2026-10-17 added the turns of Gifford's
   Fig. 6 (the presets now take the micrograph's leg axes, with turn walls a median 1.4 µm from the mask) and Xia's Fig.
   3g wall (the variant with F 0.94072 and D 6.5705 µm read off Figs. 3g and 3a, which panel f, recorded at its own
   scale, does not support along leg 7). The review also showed that the design-level Gifford Fig. 4 drawing starts its
   side channels narrower than D (LIMITATIONS 26).

7. The items below record what the papers' devices show beyond rc6, as found in rc7 (register `RC8_DISCREPANCIES.md`,
   rc7's with the RC8 entries, §§ 1 to 4 and 9; every item built and traced; the numbers from the rc7 records, and the
   Fig. 4 bullet's from the RC8 record and traces; two independent checks against App. A, App. B, and the glossary,
   `review7/foundation_check*/`).
   * For Dinh's square head corners (Fig. 2B), rc7 adds the corner term c_corner, a new pair of 2-D solves (Hele-Shaw
     and planar Stokes, `provenance/rc7_scripts/corner/`), with the solver checked against rc6's own sharp-corner fits
     and exact curved-duct solution. The inner (wall-end) rounding moves the term from −0.24 to −0.56 w, whereas the
     outer rounding does not move it measurably (to r_out/w 0.307). The term is the glossary's L_eq
     ("bend-equivalent length", an external input that this work supplies, G:1448), but it is computed rather than
     "calibrated once at low Re" (spec v5 S:409), so no validator in the record discharges it. It is also negative
     against the corner-point centerline, so spec v5's "typically 2–6 %" (S:407), a statement about U-turn folds, does
     not transfer. Its band is the head's widest, with 121 rows beyond ±1% (111 head rows and the last 10) and 4.3% of
     F, against 11 with rc6's arcs (1 head row and the last 10; LIMITATIONS 8). Every head corner is drawn rounded, at
     least at the 5 µm fillet, and the head gaps' entrance corners are drawn at 11 µm, in the drawing only (c_ent is
     external, spec v5 S:756).
   * For Xia's taper (XIA-1), the parallel-plate shear 6Q/(a b²), w_c ∝ Q_c, places leg 3's center widths on panels d
     and e within an image pixel, where the exact-duct form misses by 11 to 20 µm, and with it panel e's two
     registrations agree (the taper was chosen with those panels in view, so they are consistency rather than a test).
     It is a reading of the figures and a reproduce-mode departure from spec v5 S16-b (S:456, S:606), and it has not
     been shown to be Xia's method, because Xia's text does not state the shear coefficient (the formulas it cites are
     Yang et al.'s, whose Eqn 3 the glossary identifies as the resistance series, G:2036 to 2038). The true peak shear
     via C(α) rises 17.4% along that taper and is 1.34 to 1.57 times the parallel-plate value
     (`validation.gates.shear_peak`). Panel f is not matched (register INV-4, open). On the shipped
     `xia2016_rs7_aligned`, its outer side-channel walls lie 39.9 / 36.0 µm inside the mask's (left / right, at panel
     f's 1.774 µm/px), which is more than on rc6's variant (LIMITATIONS 24(d)).
   * Gifford's Fig. 4 (INV-2, GIF4-1; RC8) is reproduced as the 'nested reading' (`gifford2014_fig4_aligned` since RC8).
     This is the printed design with every printed value, the legs nested, Gifford's start rule at his 10 : 1 limit
     (8 µm, App. A's other start mode, with a mass defect of 0.0022; rows 1 to 3 above D, row 1 at 2.0 D, 0.5% of F; the
     8 µm row-1 channel is narrower than the 10.2 µm particles that those rows pass), and the inset's measured gap
     (17.2 µm). Four lines of evidence support the nesting. They are the caption ("the second and fourth legs"); the
     schematic's topology (its "not drawn to scale" is about scale); the inset's two arrows, which point the same way
     (the weakest line, because they may both mark leg 4); and the taper direction, which is the inset's own test of the
     layout. In a CIF leg, the side channels widen along the flow, so the sign of a leg's outer width slope is its flow
     direction, and through the trace's registration the neighbor's slope on the image has the nested leg 2's sign and
     the opposite of the path-order leg 3's. This is V13, a sign test whose standard error is the larger of a Newey-West
     one (lag 100 µm) and the independent-sample one. The nested reading gives PASS, with the image's slope
     −0.00240 ± 0.00059 µm/µm (4.1 SE) against the mask's −0.00257, at 130 of 130 eligible positions at a ±30 µm window.
     The path order gives FAIL, with −0.00254 ± 0.00062 (4.1 SE) against +0.00301, at 117 of 117. Leg 4's own array part
     gives PASS, with −0.00469 ± 0.00040 (11.6 SE) against −0.00800, at 80 of 80 at ±12 µm
     (`provenance/rc7_traces/code/taper_sign.py`). On leg 4, the claim has a median of 0.61 µm, a p90 of 1.40 µm, and a
     coverage of 0.999, which is not a blind test, because the gap is read on the same inset. State what is not
     reproduced. Leg 2's side channels beside the inset read 17.5 to 18.1 µm per side wider on the image than the
     printed single-f schedule draws them (105.1 to 106.7 µm drawn on the rows that pass the inset). The image's 122.6
     to 124.8 µm, which are the drawn widths plus the walls' measured offsets (the width profiles read 124.5 /
     125.0 µm), are what the schedule draws about 280 rows (12.7 to 13.0 mm) further on (1.17 ×). A registration or
     wall-finding error cannot be that large, and no reading with the printed single f and a start inside the beads'
     bracket explains it (open). The readings that fit leg 2's widths are rejected. Reading B (the legs in path order,
     rc7's device, with its leg 4 end fitted to those widths) is rejected by the taper direction, the caption, and the
     schematic (and by the arrows, read as one per leg). Reading C (a 51 µm start) is rejected by the inset's own beads,
     because laminar layering puts the earliest-extracted fluid against the side channels' outer walls, and on leg 4
     that band is nearly bead-free (2 beads, where any start that admits the inlet suspension puts ≥ 98), which rules
     out C's start flow (App. A's "fictitious inlet side flow", 2σ0 = 0.485). Reading D (a larger f over legs 1 and 2)
     contradicts the single printed f_gap. The sealed design start (`gifford2014_fig4_design`) puts every row at D, but
     its 4.86 µm first side channel at 80 µm depth is 16.5 : 1, against the "10 : 1 aspect ratio limitation" Gifford
     prints (p. 4503), whereas the 8 µm start keeps that limit. State the 8 µm as an inference from that printed limit
     (Gifford prints no first width and no floor). The beads bracket a straight start at up to about 12 µm and reject
     17.2 µm (the inset's own gap) and the lab's later 19 to 20 µm floor (Lezzar et al. 2022). A 20 µm floor keeps a
     sealed start only by folding the narrow rows into a head, which Gifford did not draw. This is the fold20 variant, a
     path-order figure, built at the RC8 code with exit 0, N 3,546 (a 59-row head), and every check passing, with no
     row's cutoff above the 8.3 µm beads and, on the inset, a claim median of 1.10 µm and a p90 of 2.51 µm at its
     20.5 µm gap, as the design preset's (LIMITATIONS 28).
   * For Dinh's length (INV-1), Figs. 2A and 3A put the linear section at 514 to 642 rows (23 to 29 mm). Holding D
     2.17 µm (the printed f* at 75 µm) at every gap draws 2,482 linear rows (112 mm) on the assumed constant-shear
     taper. On any taper, the linear rows after the reconstructed head need ≥ 1,764 rows (79.6 mm) at the 1-D outlet
     identity's F (0.5593 at the printed 44 / 52 µm) and ≥ 1,495 rows (67.4 mm) at the measured 0.505, and every row
     counted, with none above f_max, needs ≥ 1,897 rows (85.6 mm; section 4;
     `provenance/rc6_scripts/dinh_length_bound.py`). With the head at f* up to the printed seam, App. B's reading, the
     device needs ≥ 92.1 mm and ≥ 80.0 mm, the head counted at the linear pitch.
     - The meander cannot shorten the device by changing the row count, because N is mass balance (A:1121 to 1132, 1 − F
       = Π(1 − 2 f_i)), and resistance decides only what each side-channel segment must look like to draw its f_i
       (A:1148, A:1218 to 1219). The meander does set the head's own length, however, because a folded row needs Λ ≥ 2p
       (G:451; spec v5 G8, p = GS + max(F_min, 2R_min), 94 µm at the presets' 19 µm floor; Dinh's 27 µm walls draw
       110 µm per head row, against the 45.1 µm linear pitch), and the closure sets where folding stops (spec v5 S:374,
       G11). So, at fixed N, D, and F, a tighter fold pitch (walls at the floor would save 16 µm per folded row), a
       finer handover width (for Dinh, a 19 µm handover instead of the printed 21 µm folds 69 fewer rows and is 4.5 mm
       shorter), or a wider head shortens the device, whereas nothing brings it below N × (G + g_p). An under-resistant
       fold over-extracts locally, and the next rows compensate, with N and F unchanged (S:757). Only a design that runs
       the head above f* on purpose needs fewer head rows (about ln(1 + 2σ_seam) / (−ln(1 − 2f)), which is 279 at f* and
       75 at 3.7 f*), at the price of D in those rows, which is what Dinh's 18 measured head units show (about 3.7 f*,
       cutoffs about 4.2 µm).
     - There are two readings, the printed one (`dinh2024_v1_depicted`: N 2,701, 148.2 × 3.16 mm, array 136.1 mm with a
       25.6 mm head) and the figure-aligned two-cutoff one (`dinh2024_v1_figure_aligned`: N 1,093, 874 of them linear,
       75.7 × 3.16 mm, array 63.6 mm with the same 25.6 mm head; Dinh's chip 75.8 mm). The closure at the printed seam
       and outlet widths fixes the linear section's total extraction (A:1126 to 1132, A:1141, A:1148), and App. A
       R13(a)'s invariant lamina (A:1294 to 1303, Dinh's own rule, B:2958 to 2962) with the assumed taper turns a row
       count into a cutoff. The final rc7 (A. H. Routt, 2026-10-19) fits that cutoff to the printed length, so that at
       3.67 µm the drawn device, with the depicted head and the inlet and outlets included, is 75.7 mm (the fit moves
       about 0.09 µm per 2 mm). Dinh's beads (App. B P-C7.2) agree with it. Against Gifford's measured soft edge (his
       Fig. 5), their passage, 3.6 to 14% of the fluid's, implies D from 2.3 to 3.9 µm, whereas Gifford's curves would
       pass 29 to 50% at rc7's first reading, 4.29 µm (the figures' 642 rows), and 11 to 35% at 3.67 µm. Two caveats go
       with it. First, the printed length is kept at the price of the figures' proportion. Fig. 2A and Fig. 3A's
       photograph put the linear section at 514 to 642 rows (4.79 to 4.28 µm at F 0.5592) and its start at 46 to 58% of
       the inlet → middle-outlet distance (at most about 800 rows), whereas the preset draws 874 rows starting at 41.5%
       (Fig. 2A's head : linear proportion is 1.09, and the preset draws 0.65). No single reading of the figures gives
       the printed length, the depicted head, and their proportion together. Second, its narrowest linear rows run
       beyond κ_g's calibration. There, f_gap reaches 6.13 × 10⁻⁴ at 52 µm, which is above the calibration cap
       5.76 × 10⁻⁴ (G:1660, Gifford's array) on the 133 of its 874 linear rows narrower than 53.6 µm (21% of the linear
       section's extraction) and above spec v5's f_soft (5 × 10⁻⁴, D3: WARN) on 252 (`gates.f_beyond_calibration`). A
       prefactor 10% off there moves D by −0.05 to +0.02 µm. Do not present it as over-extraction c (A:809 to 811),
       because it is a different cutoff, row by row. It contradicts Dinh's "desired size cut-off throughout the device"
       and weakens App. B C7's wording for Dinh (B:993 to 996 and 2974 to 2975; proposed erratum D99). Under the frozen
       model, Dinh's text and figures cannot both hold at the 1-D F, and at the measured F they cannot both hold only
       with the head premise above.
   * For Gifford's Fig. 6, the 40 µm start predicts legs 3 and 5 from leg 1. This is a partial test, because the text's
     144 µm, which sets leg 7's end through F, was preferred over the caption's 140 µm with legs 3 and 5 in view
     (GIF6-2), and leg 7 is not a test. Its declared off-design entrance is rows 1 to 40 above D, row 1 at 9.9 D
     (35.9 µm), and 14.5% of F. Its early side flow, App. A's mass defect 2σ0 = 0.11 of the inlet flow (0.098 of the
     mass identity's 1 + 2σ0, G:545), is consistent with the 80 to 85% platelet recovery (Gifford's "yield", which is
     not App. B's F, B:4199 to 4201), and no particle image tests it.
   * The tools' items in rc7 are BUILD-1 (a photomask depended on what the same process built before it; fixed), the
     window rule of the trace metric (a half-window at most half the gap), the passage report (ADJ-C9, with the largest
     realized cutoff), the fluid section (ADJ-C5), the peak-shear report (S:606), and three independent code reviews and
     a convergence review before the freeze. The rc7 traces' limits are rc6's, except Dinh panel 1's V5 across (0.975,
     1.015; the square head prefers 1.011) and its null test (chamfer, the metric null ratio being 1.17), the Fig. 4
     design's p90 (3.5 µm), and the new masks' claims, all of which were set with the results in view. In rc7, 21 null
     ratios lay below the 1.3 criterion (14 on 11 claimed panels; six pass V2 on chamfer ratios of 1.05 to 1.28; the
     metric ratio is below 1 on Xia c, both masks, and on the aligned mask's Xia e), so those panels' residuals are weak
     evidence (`trace_anchors.json` null_flags; LIMITATIONS 24). RC8 (A7) counts every eligible sample. Each record
     gives the metric's coverage and its censored null ratios, a claim needs coverage ≥ 0.5, and Xia Fig. 3d (both
     masks) and 3g (the depicted mask), whose metric covers 42 to 48% of its samples with censored null ratios of 0.94
     to 1.00, are recorded rather than claimed (D5). With the censored ratios counted, 33 null ratios lie below 1.3 (16
     on claimed panels). RC8 adds V13, the taper-direction sign test (Fig. 4).

8. Release 1.0.0 adds the physics audit of every photomask, with its own ruler, and the design code's reports. The
   design code's physics auditor (`CIF_v2.py`'s P1 to P10) is ported and run on the photomask file of every build, read
   back from the DXF (`opencif/checks/auditor.py`; `docs/PHYSICS_AUDIT.md`), with three checks added. These are P0 (the
   file is the drawing), R0 (an independent ruler, `opencif/checks/ruler.py`, re-measures every ordinary straight row
   from the file and holds the generator's measurement to it within 10⁻³), and R6b (the design code's acceptance check
   of the posts, on every post). The design code's gate M14 is enforced on the drawn mask, and the design code's
   mandatory reports and acceptance validators (pressure drop, ηF, transfer bias, sensitivity map, porosity, operating
   regime, steric screen, and fabrication tolerance) are recorded on every build by its own formulas
   (`opencif/checks/design_code.py`). The evidence is the port against CIF_v2's own auditor and seeded faults
   (`tests/test_auditor.py`), each report against CIF_v2's own function (`tests/test_design_code.py`), the kernels bit
   for bit (`tests/test_kernels.py`), and the audit of every shipped photomask, in the build and again from scratch
   (`docs/VALIDATION.md`). `docs/CIF_V2_REQUIREMENTS.md` lists every requirement of `CIF_v2.py` and its status.

## 4. Claims to make carefully

* The first claim concerns the physics auditor, of which App. B says
  "the physics auditor separately checks the drawn device's inlet against its own, stricter channel limit" (P8). Since
  1.0.0, opencif runs the design code's physics auditor, ported, on the photomask file of every build, read back from
  the DXF (`validation.physics_audit`, `assessment.status.physics_audit`; `docs/PHYSICS_AUDIT.md`). Say so precisely.
  The auditor is `CIF_v2.py`'s (P1 to P10), ported, with changes of formula for two reasons. opencif's kernel pair and
  designs require two of them (P1 anchors the cubic pair, and P6b uses App. A's march-outlet identity, because CIF_v2's
  closed form exists only for its own canonical taper), and the DXF file requires the others (P9 and P10 are measured on
  the polygon). P4 compares the drawn fractions with the intended network, which is the kernel for a design without
  declared departures, because a declared departure such as Xia's loop pull-in shifts the rows downstream of it without
  decaying. Two checks are added (P0, the file is the drawing; R0, an independent ruler on the straight rows), the
  design code's own acceptance check of the posts (R6b) is applied to every post, and there is one change of role. The
  change of role concerns P8, the stricter inlet limit (Re ≤ 1, a settable default; App. B gives no value), which is
  checked on the drawn inlet and recorded under its own name rather than enforced, because at that limit every lineage
  device fails at its own flow (Dinh V1's inlet Re is 4.65 at 30 µL/min; CIF_v2 runs its lineage demo at 2 µL/min),
  whereas every opencif design meets it at 5 µL/min (inlet Re 0.38 to 0.93). The envelope that the thesis claims, the
  design code's channel gate Re ≤ 10 and the Stokes gates (M14), is enforced on the drawn mask (check M14). The
  auditor's independence has a stated boundary, because its own ruler re-measures the straight rows (most of every
  device's gaps) from the file, whereas the head units, U-turns, and outlet are taken from the generator's measurement.
  App. B's sentence that the software is "warranted by their acceptance suite and the physics auditor" then holds of
  opencif's auditor, with the difference that it reports rather than enforces the stricter limit. The erratum states
  this (D106, revised for 1.0.0, into which D95 is closed; `foundation/ERRATA_NEXT_UNFREEZE.md`) as follows.
  > A:184–187, after "the design code `CIF_v2.py` (Results 10–14)":
  > "whose kernels the targets-to-photomask generator opencif 1.0.0 (code stamp `{{STAMP}}`) reproduces bit for bit";
  > B:4375–4377: "... warranted by their acceptance suite and the physics auditor (the design code's, which opencif 1.0.0 runs on every photomask file: it enforces the design code's channel and Stokes gates on the drawn mask and reports its own stricter inlet limit) rather than by a defeasible inference."

  Do not call `CIF_v2.py` "the reference implementation", because in the frozen texts that phrase names the preprint's
  notebook, and App. A's name for it is "the design code". Cite the port's validation (`tests/test_auditor.py`, with the
  same numbers as CIF_v2's own auditor on its demonstration devices and fourteen seeded faults, which exercise every
  hard and envelope check except P2 and P6b) and the
  release's audit of every shipped photomask (`docs/VALIDATION.md`).
* The last rows and the outlet need the same care. The last rows' cutoffs rest on the outlet balance, since d(D_N/D)/d
  ln R_s,out is −26 to −175 over the 50 drawn presets and examples (LIMITATIONS 15a), which is a property of CIF's outlet closure
  that the papers' devices share. rc4 left the width change at the end of the outlet stretch out of the 1-D model, and
  2-D solves put it at 0.34 to 3.49% of R_s,out on rc4's masks (for Dinh V1, 3.49%, which would put its last 20 rows
  outside D ± 1%, 8 of them without a cutoff, and move F by −0.0072). Since rc5, the closure counts that term (a deep
  and a shallow 2-D value placed by the depth weight taken from the curved-duct solution), and the model band reports
  its range. Say
  "the last rows realise D ± 1 % in the 1-D model with the stretch-end term, within the band the record gives" rather
  than "at every gap" without qualification, and name the 3-D outlet bifurcation (η) as unmodeled. Call the band what it
  is, a sensitivity scenario of stated size (the depth weight moved, c ±30%), rather than a calibrated interval. On the
  shipped designs, it contains both 2-D limits, the drawn arcs of Dinh's flare (under 1%), and the teardrop lead's first
  turn (+8 to +14% of c, solved at three width ratios and interpolated to the nine teardrop steps), but not both limits
  at every width ratio, and the ring leads' turns and the arc after Dinh's flare were not solved (LIMITATIONS 8). Dinh's
  outlet flare (Fig. 2B, fourth panel, 14.3 to 14.6° per side over the fit intervals and the straight walls' fitted
  directions, after a straight 65.4 ± 2 µm, measured at the panel's 50 µm bar read with its anti-aliased ends; rc5 read
  67 at 46 px) is drawn by the two Dinh presets at 14.5° after 65.4 µm (rc6; 67 µm in rc5) and checked against the
  figure (`tests/test_figure_anchors.py`). The lead width that it reaches is the generator's, and the center channel's
  widening over about the last 37 µm of the panel is not drawn.
* For Dinh's length, do not write that the paper is internally inconsistent, nor that there is no conflict at the
  measured F (rc5's draft said so; both reviews of rc5 found it overstated, and independent checks of the verification
  corrected further details, among them the premise and the causes). Write the conditional, with the premises first. The
  first premise is the frozen cutoff law's cubic kernel with the lamina the same at every gap, at the plateau prefactor,
  fixed by the printed f\* = 1.04 × 10⁻⁴ at 75 µm. That is App. B's own reading of Dinh's rule
  ("Under the Xia-scaled reading, which is the one Dinh describes, the lamina is the same at every gap", B:2970 to 2972;
  Dinh adjusts each gap "to maintain the desired size cut-off throughout the device", p. 916), and under it κ drops out
  of the row count. The second is that every center width is at or above the printed 52 µm, so no row extracts more than
  f_max = 2.16 × 10⁻⁴. The third is two gaps per row. The last is that every gap row is at least the 45.1 ± 0.2 µm pitch
  measured on Fig. 2B's second, third, and fourth panels (six post-column readings, 45.05 to 45.17 µm, mean 45.11, SD
  0.04; each panel's 50 µm bar read with its anti-aliased ends, 46.91 px, to about ± 0.3%; rc5 read 46 px and
  46.0 ± 1.0 µm).
  Then (1), whatever the taper, the gap rows alone need at least 1,897 rows (85.6 mm) at the F that the presets close
  (0.5593, the 1-D outlet identity at the printed 44 / 52 µm), 74.6 mm at F 0.5105, and 73.4 mm at the measured F 0.505
  (Fig. 4, "50.5 ± 3.3%", mean ± SD; p. 916, "50.1 ± 1.6%"), against "An assembled CIF device (∼3.2 mm × 75.8 mm)" (Fig.
  3A caption, p. 918), which also "included a single inlet with an integrated filter … and two distinct outlets" (p.
  916). At the measured F, that bound leaves 2.4 mm for the inlet, filter, and outlets, whose sizes are not printed, so
  it does not show a conflict there by itself, because it lets the head's rows run at f_max. With App. B's reading of
  the head, at f* up to the printed 21 µm seam (B:3208 to 3210), the seam-anchored bound is ≥ 80.0 mm at the measured F,
  so under that one further premise the conflict holds there too (the rc7 independent foundation check; § 3 item 7
  states it with that premise). (The presets' reconstructed heads draw their first rows above f_max, a departure from
  the constant lamina, and with those rows credited the bound is 72.6 to 73.0 mm, which leaves up to 3.2 mm.) (2) The
  bound needs the center channel at 52 µm in every row, whereas Dinh narrows it "progressively" in the linear section
  (p. 916), "from 75 μm near the inlet to 52 μm near the outlet" (p. 917). Any such taper needs more rows. At F 0.505,
  it needs 105.8 mm of rows for a width linear in the row index, and 125.5 to 128.6 mm of head and rows (99.9 mm of
  linear rows) for the presets' own schedules. The rows fit the printed length only if they extract on average at least
  96.9% of f_max, as a constant 52.8 µm center would (for a linear narrowing followed by 52 µm, the narrowing within the
  first 10.2% of the rows, 11.9 to 13.6% crediting the presets' heads, and 8.8 to 11.7% over the pitch's ± 0.2 µm)
  (`provenance/rc6_scripts/dinh_taper_alternatives.py` and `dinh_schedule_at_F.py`, rc5's scripts ported to take the
  preset's own pitch, which reproduce rc5's numbers exactly at 46.0 µm).
  So, under these premises, the model conflicts with the printed footprint at the presets' F for any taper, and at the
  measured mean F for a progressive taper at the measured pitch, whereas at one SD below the mean (0.472) a linear
  narrowing confined to the first 39% of the rows would fit. It is a conflict between the premised model and the printed
  footprint rather than among the printed numbers. Its cause is not known. It may lie in the kernel's lamina-flux
  relation in Dinh's 140 µm deep channels, the constant lamina, the pitch measurement, or a taper or design rule that
  the paper does not print. The cutoff law below its plateau (Dinh's rows wider than 62 µm) is not a way out, because
  holding D itself constant with the frozen Hill curve there either makes the device longer (the every-row bound alone
  is 132.1 to 142.3 mm at the cutoffs that return the printed f\* with that curve) or does not return the printed f\*
  (at 1.805 µm, the lower end of App. B's band for a constant lamina, and at 2.1742 µm), and none of these readings
  removes the conflict at the measured pitch (a linear taper still needs 91.6 mm at the shortest). The outlet identity's
  F does not explain it, because at the measured F the bound's conflict goes, whereas the taper's does not. Fig. 2A's
  main panel cannot settle it, because read at the arm span as the 3.2 mm width, it gives 60.6 mm against the printed
  75.8 mm, so its global scale is not established (the caption's inset is "magnified").
* On scope, App. B's thesis (Leg 5) says "Rigid spheres in Stokes flow at nearly constant viscosity". State it by scale,
  as the glossary does: "The framework is therefore not globally Stokes: it is scoped to rigid, approximately spherical particles in a Newtonian, nearly constant-viscosity fluid, with creeping-flow conditions at the gap and lamina scales that govern particle selection."
  The channel flow is laminar (Re ≤ 10, gated; 3.6 to 5.6 on the strict designs), and the 1-D model needs no more,
  because in a fully developed unidirectional channel the Navier-Stokes equation "reduces exactly" to the Stokes
  equation "whatever the Reynolds number" (App. A §Setup). Do not write that the flow is Stokes everywhere.
  Deformability, cell crowding, and whole blood are declared out of scope. Platelets are discoid, and the diameter
  entering the cutoff law is orientation-dependent, so platelet-specific validity is not established. D is a
  rigid-sphere diameter. The examples named "platelet" and "wbc" are rigid-sphere designs at those sizes rather than
  validated cell separators.
* The sentence "Every gap realises D within ±1 %" is a statement about the 1-D network on the drawn geometry, at nominal
  dimensions and the nominal junction terms. Across the junction fits' own range, 3 to 5 seam rows of every sealed
  headed design move by 8.5 to 48% of D (0.1 to 0.9% of F), and since rc5 the last 6 to 58 rows of every design whose
  outlet changes width (24 of the 50) move beyond ±1%. The band's largest cutoff change on those designs, over all its
  rows, is 0.004 to 0.76 D, and its share of F, seam rows included, is 0.02 to 1.9%, and 4.3% on Dinh depicted, whose
  square-corner term adds 111 head rows. Each record names those rows (the model band).
  The statement is not a fabrication claim, and per-gap cutoffs are very sensitive to fabrication, as the records'
  design-network scenarios show (`validation.robustness_design_network`; RESULTS.md columns). A uniform ±0.5 µm width
  bias (an illustrative perturbation; no process data is cited for its size) moves 17 to 223 gaps of a device without a
  head outside D ± 1% (up to 21% of N; cutoffs 0.84 to 1.16 D, worst at a sealed start's first rows), and 82 to 404 gaps
  of a headed device (declared rows, outside at nominal, counted too; up to 23% of N, on dc_D3_F0p3 and dc_D4_F0p3;
  worst at the seam rows, which amplify element errors about 70-fold; counts over both signs). One seeded draw of 10 nm
  per-segment side-width noise (a single realization rather than a distribution) moves 0 to 2,005 gaps. Say
  "the mask realises D ± 1 % at nominal dimensions", because the fabricated chip's per-gap cutoffs will spread by
  several percent at the entrance and the seam. Gifford names the same effect in his devices
  ("imperfections ... can serve to elevate the effective value of fgap in the affected areas", p. 4500). A tolerance
  budget (the bias at which the share of F extracted within ±2% or ±5% of D falls below 95%) is future work.
* κ_g itself is known to ±3.4%, so a ±1% design tolerance is tighter than the law's calibration. κ_g was fitted on
  Gifford's Fig. 3 thresholds, "the filtration fraction per gap below which the beads of a specific diameter are consistently maintained in the central flow channel"
  (p. 4500), measured on "short-length (several centimeter long) test devices" (p. 4499), which is an onset-of-loss
  criterion. Gifford's held-out Fig. 5 devices are long (1,575 to 7,161 gaps as reproduced) and report cumulative
  enrichment. This project's reading (§2.5; `provenance/gifford_fig5_d50.py`, three interpolants × three full-retention
  ceilings) gives the retained fraction at the design cutoff R(D) = 0.39 to 0.42 (5A), 0.44 to 0.54 (5B), 0.38 to 0.48
  (5C), and 0.17 to 0.30 (5D), and d₅₀/D = 1.09 to 1.14, 0.98 to 1.03, 1.01 to 1.04, and ≥ 1.10. Read correctly, this
  says that a long device loses about half of the beads of size D, and that near-complete retention needs particles well
  above D, as Gifford writes: "when a value of fgap is chosen that generates an array with a cutoff size several microns below a given particle size, ... the vast majority of those particles are indeed retained"
  (p. 4500). It is not a test of κ_g at its ±3.4% level (as the adjudication review's E3 states, a d₅₀-within-CI test is
  ill-posed under an onset criterion), so do not write "5B within 2 % of D" or "5A three times κ_g's interval above D".
  5C and 5D are outside the envelope, and 5D's d₅₀ depends on the third digit of the enrichment ceiling. The three
  ceilings are not three readings of one number. Gifford's text gives the side-to-center flow ratio as "~23 : 1" and the
  caption a "theoretical 24.4×" enrichment, which agree within the text's rounding (24.4× is 23.4 : 1), and the script
  also carries 23 : 1 read exactly (24.0×) and the duct law at the stated widths and depth (23.5×). The 10 µm
  enrichments of 5A and 5B exceed the duct-law ceiling by 2.5 to 2.9% (R 1.025, 1.029) and the exact 23 : 1 one by 0.4
  to 0.8% (R 1.004, 1.008), and only the caption's ceiling keeps them below 1 (R 0.988, 0.992). An excess of a few
  percent is within the plausible scatter of a Coulter count, so say that the data favor the caption's ceiling rather
  than that they rule the others out.
* The Fig. 4 start width (8 µm) is an inference, so quote its source with its context. Gifford's only statement is about
  deep devices, "The relatively large (>15 μm) typical pore size for the microfluidic method we describe enables the production of correspondingly deep devices (i.e. >150 μm, assuming a 10 : 1 aspect ratio limitation)"
  (p. 4503), and applied to the 80 µm device it gives 8 µm. The same lab's later statement (Lezzar et al. 2022, Sci.
  Rep. 12:13798, Methods, "microfluidic devices with a 'minimum feature size' less than 20 µm are very difficult to mass produce … manufacturability of microfluidic devices with aspect ratio (depth:width) above~7:1 is very limited")
  gives 11.4 µm at 80 µm by its aspect ratio, inside the beads' bracket (up to about 12 µm), and 20 µm by its feature
  size, which the beads reject and the inset's own 17.2 µm gaps contradict for this device. So say that Gifford printed
  neither a first side width nor a floor, that 8 µm follows from his stated limit, and that the inset's beads allow up
  to about 12 µm.
* Fig. 4's open item is part of the reproduction rather than a footnote. State it beside the claim, with leg 4
  reproduced (median 0.61 µm), leg 2's side channels 17.5 to 18.1 µm per side wider on the image than drawn (register
  GIF4-1, open), and the layout supported by four lines of evidence, of which the taper direction is the only
  quantitative test made on the inset itself.
* F is the 1-D filtrate fraction. The glossary's η = 0.903 (a measured 3-D outlet correction rather than an external
  free parameter) is not applied; Dinh's network F is 10.7% above his measured 0.505 by construction (App. A R8; App. B
  P-C4.2).
* Most of the published reproductions sit outside the validity envelope (11 of the 18 presets exit 3). They are Xia
  RS-7 (u_min 0.021 to 0.023, M13 0.17 to 0.26, Re ≈ 61 at 0.8 mL/min, Stokes numbers up to 3.6 × 10⁻²), Dinh V1
  (u_min 0.029), Gifford Fig. 6 (u_min 0.029), and Gifford Fig. 5C and 5D on gap dominance (M13 0.083 and 0.138, as
  their records state; 5D also on the lamina Stokes gate, Re_δ 1.02 × 10⁻²). So, of the Fig. 5 devices, only 5A and
  5B are inside the envelope. Of the aligned presets (rc7's, and RC8's nested Fig. 4 reading), Gifford Fig. 4 aligned
  is inside it (u_min 0.052, M13 0.048 with the inset's 17.2 µm gap), as the Fig. 4 design is (u_min 0.052, M13 0.031,
  f ≤ 3.3 × 10⁻⁴; reported above the process aspect-ratio limit, 16.5 : 1 at its 4.86 µm first side channel). Dinh figure-aligned is inside it by its record (exit 0, because its linear section is on the plateau,
  u_min 0.049), although not in substance, because its declared head rows lie below the plateau (u 0.029, 7.2% of F) and
  its narrowest linear rows run above f_soft and κ_g's calibration cap (f up to 6.13 × 10⁻⁴; § 3 item 7), and both are
  reported rather than gated. Gifford Fig. 6 aligned (u_min 0.029) and Xia aligned (u_min 0.021) are outside, like the
  presets they align. The tool reproduces and measures them and reports a cutoff band there.
* The fluid in every shipped preset's Reynolds and Stokes gates is water (ρ 1,000 kg/m³, μ 1 mPa·s), the default of
  rc7's `fluid` section (ADJ-C5), which sets ρ and μ for those gates only (the schedule does not depend on them). App.
  B's thesis notes that the channel gate binds for water at 90 µL/min (Re ≈ 16), and that for platelet-rich plasma at
  that rate, Dinh's V1 operating point, passing turns on a viscosity that no paper states.

## 5. Notation map (code → glossary)

| code / record | glossary | note |
|---|---|---|
| `process.gap_um` | G, the link width (along the flow) | |
| `process.post_length_um` | g_p, the glossary's "post width" (the along-flow wall segment) | Λ = G + g_p |
| `process.post_width_um` | ℓ_g, the glossary's "link length" (the wall the gap crosses) | the code's names predate the glossary's, and the mapping is fixed here and in `docs/OPTIONS.md` |
| `u`, `u_min`, `plateau_u_min` | d̃ = d_crit/w_c, the reduced size ratio; u_lo = 0.035 and u_ff = 0.012 keep their names | App. B C8 itself writes u = d/w_c |
| `wc[i]` | w_c(i) | the code's `wc[i]` is the width at gap i, tapered linearly between gaps, whereas the glossary's w_c(i) is the segment after row i, which is an O(f) convention difference (App. A R13(b)) |
| `f`, `f*` | f_gap, f\* | |
| `Rc`, `Rs`, `Rg` | R_c, R_s, R_gap | |
| `S0` | σ₀ (external; the floor start's fictitious inlet side flow) | App. A R10 remark |
| `filtrate_to_retentate_ratio` | none | F/(1−F), which is not the glossary's FR |

## 6. Limitations and future work

`docs/LIMITATIONS.md` lists them with their size and the record field that reports each. None needs new computation for
the dissertation; each is stated as a limitation, and the 3-D work is future work.
