# opencif ↔ Appendix A, Appendix B and the glossary: alignment matrix

This matrix lists every statement of the frozen foundation that the generator implements, relies on, or must not
contradict, with the code that implements it, the check that pins it, and its status. We audited it against the
frozen v114 sources (2026-09-29; `foundation/`: `APPENDIX_A_revised.tex/.pdf`, `APPENDIX_B_current.tex/.pdf`,
`glossary_official.tex/.pdf`, `PREAMBLE_FINAL.tex`, and the design specification `Ladder_Array_Design_Spec_v5.md`; md5s
in `foundation/README.md`, asserted by the trace test). rc1 was audited on 2026-10-09 against v110, the previous freeze.
For rc2, the v110 → v114 differences (`foundation/v110_to_v114_*.diff`: A 20, B 175, glossary 192 changed lines) were
read in full. App. A's equation numbers are unchanged, no statement that the code relies on changed its meaning, and
AG2, AG9, and AG10 gained direct support in the v114 glossary. Every finding of the audits was fixed in code or docs or
is stated in `docs/LIMITATIONS.md`, and independent reviewers re-checked the result. The row identifiers are AA-n
(Appendix A), AB-n (Appendix B), and AG-n (glossary), so they cannot be confused with App. B's claims C1 to C9, its
script stamps "B5"/"B6", or the Prospectus items C7/C8. In the tables, "trace" = `tests/test_appendix_trace.py`, "neg" =
`tests/test_negative.py`, "units" = `tests/test_units.py`, "fd" = `tests/test_bend_fd.py`, and "Q1/Q2" = the build-time
checks on every drawn mask.

Each row has one of five statuses. A row is ALIGNED when the code implements the statement as written and cites it
correctly, and GENERALISED when the code implements a stricter or more exact form, in which case both are cited. A row
is a LIMITATION when the statement is not implemented and is stated in LIMITATIONS.md, OUT OF SCOPE when it is theory or
interpretation that the generator does not use, and N/A when it belongs to the other kernel pair or to a closure that
opencif does not use. The audit's two failure statuses, MISLABELLED (code right, citation wrong) and CONTRADICTED (code
or doc says what the frozen text rules out), no longer apply to any row, although an evaluator re-auditing a row against
`foundation/` would use them. (rc1's AB18 and AB26 were labeled GENERALISED but were CONTRADICTED, because rc1 failed a
design over its chip area, where App. B's thesis says that the design tool
"reports, without refusing, a design over its length budget (C8, CQ5)" and CQ5 describes that report as made
"without failing". rc2 reports it. rc1's AG7 named γ_GS and f_soft without a disposition, and rc2 ports both. These are
F2 and F3 of the review of 2026-10-10.)

## 1. Appendix A (equation numbers as printed)

| # | source | statement | code | check | status |
|---|---|---|---|---|---|
| AA1 | §Setup, eqs. 1 and 3 to 5 | steady, laminar, unidirectional flow, for which the Navier-Stokes equation "reduces exactly" to the Stokes equation "whatever the Reynolds number, as long as the flow is laminar and fully developed"; creeping flow at the gap; f = flow fraction of the near-wall lamina; x = δ/w_c | `CIF_v2.g_cubic`; the 1-D network's resistances do not depend on the flow | none | ALIGNED |
| AA2 | R1, eq. 10 | f = 3x² − 2x³ | `schedule.f_of` → `kernels.g_cubic` (= `CIF_v2.g_cubic`, bit for bit) | trace "A R1 (eq. 10)" | ALIGNED |
| AA3 | R2, eqs. 18 and 20 | the leading law δ ≈ w_c√(f/3) and the exact root | `schedule.x_of` → `kernels.x_of_f(..., 'cubic')` (= `CIF_v2.x_of_f`, bit for bit) | trace (round trip, leading order, and the implementer roots x(10⁻⁴), x(5.8 × 10⁻⁴), and x(0.1)) | ALIGNED |
| AA4 | R3 and R4 | the CIF and DLD cubics unified | none | none | OUT OF SCOPE |
| AA5 | R5, eqs. 34 and 35 | κ_steric = 2 and the κ_app decomposition | absorbed in κ_g | none | OUT OF SCOPE |
| AA6 | R6, eq. 36 | the operating law d_crit = κ w_c √(f/3) (leading order) | δ\* = D/κ_g, f(i) = g(δ\*/w_c) on the exact root, as R13(a) does | units; Q2 on every mask; trace "exact root vs leading order ≤ 0.47 %" | GENERALISED (exact root) |
| AA7 | R7, eqs. 41 and 42 | the corrected Xia rule, whose exact form holds δ\* | `schedule.build` | Q2 on every taper; trace "lamina at every row is δ\*" | ALIGNED |
| AA8 | R8, eqs. 44 and 45, P1 to P3 | F = 2R_c/(R_s + 2R_c) under outlet pressure equilibration | not used, replaced by the exact outlet closure (`schedule.build`/`reclose`) | trace "exact outlet vs A R8 eq. 44: ≤ 4 f_gap"; "network delivers the design F" | GENERALISED |
| AA9 | R8, eq. 46 | ∂F/∂f_gap = 0 at fixed outlet geometry | the stop rule on F with the exact closure (the outlet is sized per design) | none (input independence is AA24) | N/A |
| AA10 | R8, eq. 47 | R_s/R_c = 2(1−F)/F | not used (outlet from the closure) | none | N/A (DISSERTATION_NOTES §2.2) |
| AA11 | R8, last paragraph | finite networks hold P3 only approximately | the exact closure does not assume P3 | Q2 | ALIGNED (cited as context) |
| AA12 | R9, eqs. 48, 53, and 55 | K_depth, F_lam, and c₂ | not used (cubic pair) | none | N/A |
| AA13 | R9 corollary, eq. 59 | one kernel pair at all three call sites | the cubic + κ_g at `x_of`, `f_of`, and δ\* = D/κ_g; `design.kernel_pair` in every record | trace "pairing rule" and "kernel declaration" | ALIGNED |
| AA14 | R10, eqs. 60 to 67 | the march and F = 1 − P(N) | `schedule.build` | trace "implementers march example" (P 0.367843, σ 0.859277, F 0.632157) | ALIGNED |
| AA15 | R10 | the mass identity P + 2S = 1 + 2σ₀ | `schedule.build` | trace (asserted < 10⁻¹², observed 4 × 10⁻¹⁴) | ALIGNED |
| AA16 | R10 σ₀ remark; eq. 67 | F = 1 − P(N) for any σ₀, where a floor start's σ₀ is a fictitious inlet side flow | a floor start stops on the F that its drawing delivers, 1 − P(N)/(1 + 2σ₀) (= eq. 67 at σ₀ = 0) | units "floor start stops on the DELIVERED F" | GENERALISED |
| AA17 | R10, eq. 69 (with eq. 68, Gifford eqn 4, eq. 64) | the PE closure σ(i) = R_c/R_s | the exact gap-drop closure (spec v5 §13, equivalence-chain identity), which reduces to eq. 69 | trace "gap-drop closure: network gives every designed f" (≤ 5.7 × 10⁻¹⁰) | GENERALISED |
| AA18 | R10, eq. 70 | the inverse map recovers f(i) | replaced by an independent network solve (`measure.solve3`) | trace; Q1 on every mask | GENERALISED |
| AA19 | R11, eq. 75 | 1 − P(N) = 2σ(N)/(1 + 2σ(N)) | the march satisfies it | trace "A R11 (eq. 75) to rounding" (8 × 10⁻¹⁵) | ALIGNED |
| AA20 | R11, eq. 76; implementers step 7 | 1 − P(N) = 2R_c/(R_s(N) + 2R_c) with the PE outlet | not met by design, because the exact outlet differs by the outlet row's gap drop | trace "exact outlet vs A R8 eq. 44: ≤ 4 f_gap" (3.3 f_gap) | GENERALISED |
| AA21 | R12, eq. 80 | 2φ_crit = F | none | trace "transport lemma" (`CIF_v2.critical_phi`, ≤ 1.4 × 10⁻¹⁴) | ALIGNED |
| AA22 | R13(a), eq. 81 | the invariant-lamina schedule | `schedule.build` | trace "lamina at every row is δ\*"; Q2 | ALIGNED |
| AA23 | R13(b), eqs. 82 and 84 | N·D² invariant at fixed F | `schedule.plan` | trace (a 0.51% spread over D from 3 to 7 µm, with D 3 just below the plateau, and 0.92% from eq. 82) | ALIGNED |
| AA24 | R13(c), eq. 86 | the lower-triangular Jacobian, by which D fixes δ\* and F fixes N (on the plateau) | `spec` → `schedule.plan` | trace "changing D changes N not F", "changing F changes N not the cutoff" | ALIGNED |
| AA25 | R13(b), Ψ taper | the canonical taper | not offered | none | LIMITATION 16 |
| AA26 | R14, eq. 88 | the rectangular-duct resistance | `schedule.rr` → `kernels.r` (= `CIF_v2.r`, bit for bit) | trace (401-term series, 1.9 × 10⁻⁹; square duct 28.454154) | ALIGNED |
| AA27 | R15, eqs. 91 to 93 | the PE imbalance decay | not needed (exact closure); M13 guards | none | N/A |
| AA28 | "For implementers", step 3 | stop if δ\*/w_c ≥ ½ or f ∉ (0, ½) | `spec.to_engine` SpecError; `schedule.f_of` ValueError; an unplannable closure is a SpecError | trace; neg (D 250 µm and D 100 µm on 100 µm) | ALIGNED |
| AA29 | "For implementers", test values | the roots, square duct, and march (the code's kernels); c₂, F_lam, the series coefficients, and the ladder decay | as above | trace (the ones that the code uses) | ALIGNED / N/A |

## 2. Appendix B

| # | source | statement | code | check | status |
|---|---|---|---|---|---|
| AB1 | §Notation | κ_app is constant only on the plateau u ≳ 0.035 | `build.gates` `plateau_ok`; strict refusal | neg "u_min 0.032 → REFUSED" | ALIGNED |
| AB2 | C1 | the cutoff set by the lamina, and κ_steric = 2 | embedded in κ_g | none | OUT OF SCOPE |
| AB3 | C2 | δ = w_c√(f/3) (leading term), with κ fitted on the parabolic δ | the exact cubic root (≤ 0.47% from the leading term up to the cap), with the depth factor never applied | trace | GENERALISED (exact root) |
| AB4 | C3a and C3b | κ = 4.91 [4.75, 5.08] (±3.4%) | `KAPPA_G` and `KAPPA_CI` from CIF_v2 Pairing | trace (macros `\kappVal`, `\kappCIlo/hi`) | ALIGNED |
| AB5 | C3c, C7 P-C7.1, CQ3, CQ4 | below the knee, the law is implicit; the Hill fit (5.568, 0.01651, 1.791); the Dinh band [1.25, 2.17], local κ 4.08, and collapse 1.80; plateau-anchored and of low severity | `build.gates` `cutoff_band_um`, with the Hill κ at u = κ_g δ\*/w_c and the band [κ(u)δ\*, κ_g δ\*] (App. B's band for a lamina held at δ\*), with the CQ3 caveat, and none below the foot | the trace reproduces 1.25, 2.17, 4.08, and 1.80, and the gates give 1.80 at Dinh V1 | ALIGNED |
| AB6 | C3d, C6 | the κ decomposition, with β = ½ and γ = 1 | none | none | OUT OF SCOPE |
| AB7 | C4, P-C4.2 | F_closed V1/V3 = 0.5592 / 0.7189 | the Dinh V1 preset's F; `devspecs.Fcl` (V2's 0.6169 marked unverified and unused) | trace `\FclosedVone`, `\FclosedVthree` | ALIGNED |
| AB8 | C4, P-C4.3 | the network F equals the identity to 0.12% | the exact outlet closure | Q2 | ALIGNED (cited as context) |
| AB9 | C4 CQ1 | a 3-D gap-mouth "sibling" effect on f_eff (about 5% on √c), with a common 3-D origin of the residuals open | not modeled | none | LIMITATION 9 |
| AB10 | C5 | the corrected Xia scaling | `schedule.build`, exact form | none | ALIGNED |
| AB11 | C9 | the transport lemma | none | trace | ALIGNED |
| AB12 | App. B's "B5", "B6" (regression stamps of its script) | none | opencif uses neither label as App. B's stamp (its internal patch is "RP-B5"); since 1.0.0 `B6_F_real` in `validation.design_code_reports` is the design code's own label for its report of F_real = ηF, a different B6 | grep finds App. B's stamps nowhere in code or tests | ALIGNED |
| AB13 | C8 claim; thesis Leg 4 | forward substitution realizes both targets "on a device built to the resulting schedule", inside the envelope, and "This appendix reports no device fabricated to such a schedule" (v114) | `spec` → `plan` → `build` → Q1/Q2 | Q1 and Q2 on every mask | ALIGNED (a 1-D network statement on the drawn mask rather than fabrication; LIMITATION 22) |
| AB14 | C8 kernel declaration | declare the pair | `design.kernel_pair` in build, refusal, and check records | trace | ALIGNED |
| AB15 | C8 P-C8.1 | δ(i) = δ\* at every row | `schedule.build` | trace; Q2 | ALIGNED |
| AB16 | C8 P-C8.2 | N·D² is invariant (App. B's 0.29% is the F_lam pair's) | `schedule.plan` | trace (cubic pair, 0.51%) | ALIGNED |
| AB17 | C8 P-C8.3 | the decoupling over 1.73 decades of f_gap | none | `examples/decoupling/` (D 3 to 5 × F 0.3 to 0.7) | OUT OF SCOPE (illustrated) |
| AB18 | C8 CQ5 (with the thesis sentence that cites it) | N must stay free; the design code's budget check "reports the result without failing", and "`design` itself refuses nothing on length"; the thesis says that the tool "reports, without refusing, a design over its length budget" | `fit_layout` varies legs and head width to fit the chip and never N; a design over its chip area is drawn and reported ("fits the chip area", kind `budget`; the verdict ends "OVER THE CHIP AREA ... reported, not refused or failed"; `assessment.over_chip_area`), with the exit code that it would otherwise have, and arrays likewise; without a chip area, nothing is reported | neg "control: a design over its chip area is drawn and reported (exit 0, ...)" | GENERALISED (rc3 label; rc2 said ALIGNED, and rc1 failed it). App. B's budget is spec v5 G9's 1-D leg length, whereas opencif reports the 2-D chip footprint, the predicate that App. B's CQ5 motivates ("to fit its overall footprint onto a standard substrate mould", Xia), under the same report-not-fail policy; `api.Result.over_chip_area` and `assessment.status.chip_area` keep it apart from the exit code. rc4 also reports G9's own predicate, each leg's path length against the chip's leg length (`validation.length_budget`, `status.length_budget_G9`), as App. B describes the design code doing ("a design's report prints the length beside the budget") |
| AB19 | C8 CQ6 | run-outs must match ρ(N) | the exact outlet closure, with the outlet leads measured | Q2 F | ALIGNED |
| AB20 | C8 CQ7 | the plateau guard refuses Dinh (0.029, and his other designs at 0.030 and 0.032) and Xia RS-7 ("runs u = 0.023 at its 300 µm wide end ... and the guard refuses it on the same test") on the cubic pair | `build.gates`; strict refusal (exit 5) | neg (the plateau, channel Re, and lamina Stokes refusals); the Dinh and Xia records (exit 3, reproduce mode) | ALIGNED |
| AB21 | thesis | D and F are the design targets, and f_gap is internal | `spec.SCHEMA` (no f_gap target); `api.generate(D, F)` | trace "f_gap is refused" | ALIGNED |
| AB22 | thesis summary | yield = F, never particle recovery | the docs say "filtrate fraction F" | none | ALIGNED |
| AB23 | thesis Leg 5, Scope | "Rigid spheres in Stokes flow at nearly constant viscosity"; deformability, crowding, and whole blood out of scope; platelets discoid, with an orientation-dependent diameter | README, DISSERTATION_NOTES §4, and LIMITATIONS 1, worded by scale as the glossary does ("not globally Stokes", with laminar channel flow, Re ≤ 10 gated, and creeping flow at the gap and lamina scales) and as App. A §Setup explains (AA1) | none | ALIGNED (stated in the glossary's scale wording, with which App. B's own scope prose agrees, "the channel-scale Reynolds number ... is ≈ 5.2", and the gap's ≈ 5.2 × 10⁻³) |
| AB24 | thesis Scope, "the Stokes regime is gated in code" | Re_stokes = 10⁻² on the gap and lamina and Re_max = 10 on the channel, and the design code refuses either (its M14) | `build.gates` (CIF_v2 M14 at every row, its row-1 term included); strict refusal | neg "lamina Re 1.05·10⁻² → REFUSED", "channel Re 37 → REFUSED"; trace "CIF_v2 M14 literal (row 1, full inlet flow) = build.gates" at 25, 28.5, and 30 µL/min (≤ 5 × 10⁻¹⁶) | ALIGNED |
| AB25 | thesis Scope | the channel gate can bind (water at 90 µL/min), and viscosity matters for plasma | in rc7, the `fluid` section (ρ, μ; Newtonian only) sets the Reynolds and Stokes gates, water by default (ADJ-C5), and the schedule does not depend on it | `validation.gates.fluid` | LIMITATION 2 |
| AB26 | thesis (introduction and Claim) | "it refuses a design that leaves the plateau or the Reynolds gates, and it reports, without refusing, a design over its length budget (C8, CQ5)" | in the strict default, refusals (exit 5) on the plateau and foot and on the M14 channel and Stokes gates; a design over its chip area or its G9 leg-length budget is drawn and reported, and never refused or failed (AB18) | neg | GENERALISED (extra refusals, which are not App. B's. In strict mode, they come from spec v5 and CIF_v2 `design()` and validators, and they are M13 gap dominance (§18 and §19), E11 thin lamina on every row, M16 taper rate, the γ_GS pre-flight (§14.4, §18), and w_min = max(F_min, h/AR_max) on G, a drawn GS, every w_c, the terminal w_s, and, since rc4, every straight side channel (G7; G11's fold is not drawn). In every mode, they are a head that cannot seal gap 1 under `seal: require` and, since rc4, an outlet that cannot close the design at any leg count tried) |
| AB27 | thesis scope note | Prospectus C7/C8 carry no argument nodes, and they are "warranted by their acceptance suite and the physics auditor rather than by a defeasible inference" (lines 4376 and 4377; the item opens at line 4373); "the physics auditor separately checks the drawn device's inlet against its own, stricter channel limit" (lines 4349 and 4350; provenance "physics auditor (P8)") | opencif implements Prospectus C7/C8 as a derivative of the design code `CIF_v2.py`, and its acceptance suite is `tests/`. Since 1.0.0, it runs CIF_v2's physics auditor, ported (`opencif/checks/auditor.py`), on the photomask file of every build. P0 (the file is the drawing), R0 (its own ruler on the straight rows), P1 to P6b, P9, R6b, and P10 are pass conditions in every mode, and the design code's M14 is enforced on the drawn mask (check M14, envelope). P8's two criteria (CIF_v2 `_PHYSICS_SRC`, with inlet Re ≤ 1.0 at the device's flow, default 5 µL/min, and 0.05 Re D_h ≤ Λ/2) are checked on the drawn inlet at the design flow and at 5 µL/min and reported under their own name (`validation.physics_audit`, check P8) rather than enforced, because at their own flows every shipped design exceeds the Re limit (3.6 to 61), whereas at 5 µL/min all pass. The warrant sentence (lines 4375 to 4377) refers back to the auditor of lines 4349 and 4350, so it holds of opencif's auditor, which enforces M14 and reports P8. D95 (up to RC8, the auditor not run on the mask) is closed into D106, revised for 1.0.0, which names opencif and says which gate its auditor enforces (`docs/DISSERTATION_NOTES.md` § 4); the frozen text is not edited | trace "P8 (CIF_v2 _PHYSICS_SRC) ...", "P8 inlet Re in the record = CIF_v2 reynolds"; `tests/test_auditor.py` (the port against CIF_v2's own auditor; seeded faults) | ALIGNED (1.0.0; P8 reported, LIMITATION 6); D106 (proposed, revised for 1.0.0) |
| AB28 | none | the curved-duct resistance | cited as new work (the counterpart of App. A R14) rather than as an App. B addendum | units (limits); fd (finite differences, ≤ 3 × 10⁻⁷) | ALIGNED |

## 3. Glossary (cited by subsection name)

| # | source | statement | code | status |
|---|---|---|---|---|
| AG1 | (AB12) | no "App. B (B5, B6)" citations | relabeled everywhere | ALIGNED |
| AG2 | Cutoff-Law Prefactor Decomposition, pairing rule; v114, "its `design` function defaults to the cubic pair" | the cubic with κ_g; never cross | the cubic + κ_g everywhere; κ_F = 5.798 never used | ALIGNED |
| AG3 | the reduced size ratio d̃ = d_crit/w_c (moved off the letter u) | notation | the records keep `u_min`, mapped to d̃ in DISSERTATION_NOTES §5 | ALIGNED (mapped) |
| AG4 | Validity-Envelope and External Inputs, u_lo and u_ff external | provenance class | PRESET_PROVENANCE codes the plateau edge A (external) | ALIGNED |
| AG5 | Geometric Variables, with G the link width, g_p the "post width" (along flow), and ℓ_g the "link length" | symbol names | `post_length_um` = g_p, `post_width_um` = ℓ_g (mapping in OPTIONS.md and DISSERTATION_NOTES §5) | ALIGNED (mapped) |
| AG6 | w_c(i) = the segment after row i | the index convention | `wc[i]` at gap i (an O(f) difference) | ALIGNED (stated) |
| AG7a | Validity-Envelope and External Inputs, c_ent (external) | supplied | `process.gap_entrance_c`, default 1; one scenario of `validation.robustness_design_network` varies the head mouths' R_gap ±20% | LIMITATION 9 |
| AG7b | L_eq (external) | supplied | for arcs, superseded by the exact curved-duct solution (no equivalent length); for rc7's square quarter turns (`meander.turns: square`), supplied as the computed corner term c_corner (G:1448 EXTERNAL; spec v5 S:404; a 2-D computation that is not calibrated, S:409; negative against the corner-point centerline, D105; its band in the record, LIMITATION 8) | GENERALISED (new work, AB28); rc7 SUPPLIED (computed) |
| AG7c | AR_max (external; spec v5 §18 w_min = max(F_min, h/AR_max), G7 "widths ≥ w_min", G11 fold iff w_s(i) < w_min; CIF_v2 `Process.AR_max` 10) | supplied | `process.aspect_ratio_max` (default 10); CIF_v2's w_min predicates as gates checked before drawing, with G and GS (when a head is drawn) before planning and every w_c and (rc4) every straight w_s once planned (CIF_v2 `design()` and its validator G7; spec v5 G7; strict refuses; neg "G 20.5 um below w_min ... refused", "a 5.35 um straight side channel ... refused in strict mode before drawing (G7)"); on the drawn mask, the narrowest CIF channel's depth/width is a pass condition in strict mode and is reported in the others (neg control "the same in evaluate mode is drawn and reported"); G11's fold geometry is not implemented, and where it would fold, strict mode refuses | LIMITATION 7 (rc4, with G7 gated in strict mode and G11's fold replaced by the refusal) |
| AG7d | R_min (external) | supplied | `process.min_radius_um`; the corner-radius check on every drawing | ALIGNED |
| AG7e | f_soft (external; spec v5 §18 "if max(f) > f_soft: WARN" and §19 validator D3; CIF_v2 `Pairing.f_soft` 5 × 10⁻⁴) | a warning rather than a gate | `build.gates` `f_soft` and `f_above_soft`; a note in the record above it (beside the calibration-cap note, 5.76 × 10⁻⁴) | ALIGNED (rc2) |
| AG7f | γ_GS (external; spec v5 §14.4 and §18 pre-flight "ASSERT GS/G in gamma_GS"; CIF_v2 `design()` GateError; `Pairing.gamma_GS` [1.4, 1.5]) | a pre-flight gate | `pipeline.gamma_gs_check` when a meander head is drawn, where strict refuses (exit 5) and evaluate/reproduce report it, with a narrow head's narrow units noted (neg "seeded: meander GS/G 1.58 outside gamma_GS -> refused"; trace "CIF_v2 Pairing") | PARTLY (rc2; the trace test labels it "partly gated"), an applicability difference. CIF_v2 checks GS/G on every design, because its narrow side channels are always folds of width GS, whereas opencif draws no GS without a head and checks it only when a head is drawn (neg control "without a meander head ... not applied") |
| AG7g | φ_cal (external, undischargeable) | supplied | not used (absorbed in κ_g) | OUT OF SCOPE |
| AG7h | σ₀ (external) | supplied | the floor start's σ₀ (AA16) | ALIGNED |
| AG7i | η (measured, 0.903) | a realized-yield factor | not applied, although since 1.0.0 ηF is reported beside F on every record (`validation.design_code_reports.B6_F_real`), as the design code reports it | LIMITATION 5 |
| AG8 | three conditions (above-knee, run-out match, and free row count) | envelope | = AB18 to AB20 | ALIGNED |
| AG9 | entrance length, "Report L_e/Λ; do not gate on it.  Gate instead on the perturbation, which the design specification identifies with the design code's gap-dominance gate (M13)" (v114); spec v5's O(f²) error order "needs the first-order term to cancel, which this record does not show, so it is not used here as a bound" | as quoted under source | `Le_over_Lambda_L` (the Stokes-floor correlation, CIF_v2 `entrance_length`) reported rather than gated; M13 gates the perturbation (R_gap f/R_c); opencif's docs do not use spec v5's "benign ... O(f²)" (§16-c, M15); P8's naive 0.05 Re D_h is reported rather than gated (AB27) | ALIGNED |
| AG10 | Reynolds number, where the argument record's scope gate uses the flux product f·Re_channel; v114, "The design code's Stokes gate (M14) instead computes the gap's own Reynolds number, over the gap's G × h section ..., which exceeds this product by (w_c + h)/(G + h) ..., and a lamina-shear Re_δ, and in a cutoff design refuses either above 10⁻²" | as quoted under source | the Stokes gates with CIF_v2's `reynolds` in the gap section, at every row, where (w_c + h)/(G + h) = 1.2 to 2.8 × the flux product on the presets | ALIGNED (the design code's reading, the stricter one) |
| AG11 | Cross-Lineage Validation, FR = Inglis and Herman's filtration ratio | symbol | the code's `filtrate_to_retentate_ratio`, documented as not FR | ALIGNED |
| AG12 | the calibration cap 5.76 × 10⁻⁴ (`\fgapDesignCap`) | κ_g's range | `f_design_max` reported, with a note above the cap | ALIGNED |
| AG13 | F_outlet vs F_recursion | two identities | opencif's F is F_recursion delivered by an exact outlet | ALIGNED |
| AG14 | α = h/w_ch | notation | not exposed | ALIGNED |

## 4. Glossary numeric macros the code uses (all pinned by trace)

The trace test pins `\kappVal` 4.9138, `\kappCIlo`/`\kappCIhi` 4.75/5.08, `\plateauUwc` 0.035, `\footUwc` 0.012,
`\squareDuctR` 28.4542, `\dcritGGIn` 2.174, `\FclosedVone` 0.5592, `\FclosedVthree` 0.7189, `\hillKmax` 5.568,
`\hillUhalf` 0.01651, `\hillN` 1.791, `\dinhBandLo`/`\dinhBandHi` 1.25/2.17, `\dinhKappaLocal` 4.08, `\dinhCollapseLo`
1.80, `\kappOverRootThree` 2.8370, and `\fgapDesignCap` 5.76 × 10⁻⁴. Two macros are consistent but not pinned (see
RESULTS.md, envelope column). They are `\ReWorstChannel` 5.2 (opencif's channel Re on the strict designs is 3.6 to 5.6
at 25 to 30 µL/min) and `\ReGapScale` 5.2 × 10⁻³ (opencif's largest Stokes number on a strict design is 8.8 × 10⁻³,
`gifford_physics_linear`). `\kFval` and `\ctwobar` (the other pair) are not used. `\etaVal` 0.903 = CIF_v2's
`Pairing.eta`, which 1.0.0 uses to report ηF and never applies (LIMITATION 5; the trace test asserts the value against
CIF_v2's).

## 5. CIF_v2's supplied constants and tolerances

Every field of `CIF_v2.Pairing`, `Tolerances`, and `Process` is accounted for, in one of five ways that the trace test
keeps apart ("CIF_v2 Pairing / Tolerances / Process: every field accounted for"). The test asserts the field sets
exactly (9, 11, and 19 fields), so a field added, dropped, or left without a disposition fails it. The counts are
Pairing 2 gated, 1 partly gated, 3 reported, 3 consumed, and 0 other; Tolerances 5, 0, 0, 0, and 6; and Process 1, 0, 0,
0, and 18.

The gated fields are evaluated on every design, strict mode refuses on them, and each value is asserted equal to
CIF_v2's. They are u_lo and u_ff (the plateau and its foot), Re_max, Re_stokes, and gap_dom (M14, M13), eps_lam and
eps_taper (E11, M16; rc3; E11 on every planned row since rc4), and AR_max. For AR_max (rc4), the gate is w_min on G and
GS before planning, on w_c and every straight w_s once planned, and on the drawn mask, where strict mode fails (exit 2)
rather than refuses, and G11's fold is replaced by the refusal (LIMITATION 7). The partly gated field is γ_GS, a
pre-flight gate when a head is drawn, whereas CIF_v2 tests it on every design (AG7f). The reported fields are evaluated
on every design and are never a pass condition. They are f_soft and, since 1.0.0, η (ηF, the design code's F_real) and
β_cal (the transfer bias that it implies, the design code's mandatory report, which on the cubic pair does not enter N),
each with its value asserted equal to CIF_v2's. The consumed fields, which hold values that the law uses, are κ_g, its
interval, and the cubic kernel. The remaining fields are a spec input, replaced, or a named limitation. They are the
process and geometry inputs; μ and ρ_f (LIMITATION 2, water by default; rc7's `fluid` section); L_eq (the exact bend
solution for arcs, and rc7's computed c_corner for square turns, AG7b); the fold stroke l_stroke; the leg length (the
chip area, reported); m_budget (the exact outlet closure); and CIF_v2's solver tolerances eps_r, eps_inv, eps_cub,
eps_id, eps_geo, and eps_net (opencif's own checks).

The test checks each gated, partly gated, reported, and consumed value against CIF_v2's and each named LIMITATIONS item
against its text, whereas the other dispositions are stated rather than tested. A field's disposition says what opencif
does with it. That a gate is implemented is shown by the negative suite and by the trace suite (the w_c floor). The
negative suite refuses a seeded defect for every gate except the foot and the Re_gap-only refusals, which are not
seeded; M13, which is seeded only together with other gates and refused alone in the record of
`wbc_xia_process_refused`; and M16, which has not been made to fail alone on a constructible design and is seeded alone
on the gate function and together with M13 as a design. rc4 adds E11 on a planned row alone, G7 on a straight side
channel, and the outlet refusal. The gates' parity with CIF_v2 is shown by the trace suite (E11 and M16 against
`CIF_v2.march` on matched requests, just below and just above each limit; w_min against `Process.w_min`; M14 against
CIF_v2's literal row-1 formulas). (In 1.0.0, the trace test's disposition text for L_eq names both of its replacements,
the exact bend solution for arcs and c_corner for square turns, whereas RC8 left the text as it was, to keep its tests
hash.)

## 6. Design specification v5 (rc7 additions)

| row | spec v5 | statement | opencif | check | status |
|---|---|---|---|---|---|
| S1 | S16-b (S:449, S:456), S:606 | a constant-shear taper holds the maximum wall shear via C(α), "NEVER 6"; report max γ̇ via C(α) | `centre_channel.shear_formula: exact` (default, C(α)); `parallel_plates` (6Q/(a b²), w_c ∝ Q_c) in reproduce mode only (register XIA-1, `xia2016_rs7_aligned`); the peak shear via C(α) reported on every record | `validation.gates.shear_peak`; `tests/test_rc7.py` | ALIGNED (exact); DEPARTURE, reproduce only (parallel_plates; LIMITATION 16). Xia's own text names rectangular-duct formulas, saying that the taper "is determined numerically (using formulas for rectangular channels based on Yang et al.) in order to maintain a relatively constant amount of (maximum) shear experienced along each leg of the device" (Xia 2016, p. 4). It does not state the shear coefficient, so it neither states nor rules out the parallel-plate form, which is what Fig. 3's panels d and e read (register XIA-1; Review 1 of rc7, Q6) |
| S2 | S:404, S:407, S:409 (§14.6) | R_s^eff = r(GS)[L̂ + n_turns L_eq], with L_eq "calibrated once at low Re", typically 2 to 6% for U-turn folds | for square quarter turns, c_corner per corner from 2-D solves (Hele-Shaw and planar Stokes, the depth weight), which is not calibrated and is negative against the corner-point centerline | `model_scenarios.bend_model.square_corner` (its band) | GENERALISED (computed rather than calibrated; LIMITATION 8) |
