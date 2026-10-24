# CIF_v2.py requirements and their status in opencif 1.0.0

Read-only audit, 2026-10-22.  Sources: `reference/CIF_v2.py` v2.1.0 (8058 lines, md5 0e0961a5..., asserted by tests/test_kernels.py:36), its embedded `_PHYSICS_SRC` (ladder_physics.py) and `_MEANDER_SRC` (meander_pattern.py), read in full; the opencif 1.0.0 working tree (`opencif/__init__.py` `__version__ = 1.0.0`).  CIF_v2 was also imported read-only to list its registry: `VALIDATORS` holds 162 entries.  opencif line numbers are in the 1.0.0 working tree.

## How CIF_v2 checks things (this decides what "port" can mean)

* **On the design being built** CIF_v2 enforces only section A: `design()` and `march()` raise on the plateau, gamma_GS, w_min, the lamina and taper limits, M13 and M14; the constructors and mask builders raise on bowtie posts, impossible U-turns, routability and (v2) composition interference; the unified compiler adds four pre-flight gates (section E, CU1-CU4).  `design()` also emits mandatory reports (section B).
* **`validate()`** (162 validators, section D) is a self-test of the implementation on CIF_v2's own fixtures (`_Ctx`, `reference_request()`, demo specs A/B/C).  It never sees a user's design.  Its own printout says so: "This proves the IMPLEMENTATION is correct.  It cannot prove the MODEL is."
* **`run_acceptance()`** passes when `goldens and G14 and jsuite` (8041-8043) hold on CIF_v2's optimizer devices; the D-suite, `V_operating_regime`, `V_steric` and `V_fab_tolerance` are written to its JSON but are not in `ok`.  `V_outlet_divider`, `J5`, `J10` and `self_test()` are defined and never called.
* **The physics auditor** (`audit()`, P1-P10 with P6b) runs on CIF_v2's meander demo patterns; Z1, Z2, T1, W3, W4, Q1, Q5 and R2a require every P check to pass, and `compile_unified` refuses when P3 or P4 fails.

Statuses: **PORTED-GATED** (the same check is a pass condition or refusal in opencif), **PORTED-REPORTED** (the same check is computed and recorded, never gated), **EQUIVALENT** (same purpose, different method: the evidence says how), **NOT-APPLICABLE** (concerns geometry or machinery opencif does not have, or CIF_v2's own fixtures), **MISSING**.  "level": where opencif's counterpart lives: build (every design), test (tests/), both, or -.  "?" marks a status I am not certain of (listed again at the end).

CIF_v2 roles: REFUSE (raises on the request), TRIP (AssertionError tripwire), REPORT (Device.reports/report()), VAL (validate() pass condition on CIF_v2's fixtures), PRE (compile_unified pre-flight refusal), AUDIT (audit() check), ACC-PASS (in run_acceptance's ok), ACC-REP (acceptance JSON only), DEF (defined, never called), SELF (embedded self-test).

## Counts

| section | rows | PORTED-GATED | PORTED-REPORTED | EQUIVALENT | NOT-APPLICABLE | MISSING |
|---|---|---|---|---|---|---|
| A | 35 | 20 | 0 | 8 | 7 | 0 |
| B | 14 | 0 | 6 | 1 | 2 | 5 |
| D | 162 | 33 | 5 | 63 | 59 | 2 |
| E | 46 | 10 | 4 | 23 | 5 | 4 |
| F | 23 | 8 | 2 | 9 | 4 | 0 |
| G | 2 | 0 | 0 | 0 | 2 | 0 |
| **all** | **282** | **71** | **17** | **104** | **79** | **11** |

Section A+B are what CIF_v2 enforces or reports on a design: 49 rows.  The 47 supplied constants (section C) are tallied separately below, since most of them are the thresholds of the rows above.

## A. Gates, refusals and tripwires applied to the design being built (design(), march(), constructors, mask builder)

These are the only CIF_v2 checks that act on a user's request in `design()` / `build_mask()`: everything in section D runs on CIF_v2's own fixtures.

| ID / name | CIF_v2 line(s) | what it checks, threshold | role in CIF_v2 | opencif 1.0.0 status | level | evidence (file:line, record field) |
|---|---|---|---|---|---|---|
| A1 plateau guard (S12) | 1222-1226 | max over 401 t of W(t) <= D_target/u_lo, u_lo = 0.035 (skipped for a prescribed f schedule) | REFUSE (WindowError) | **PORTED-GATED** | build | build.py:864-865,960 `gates.plateau_ok` (u_min = D/max w_c >= 0.035); pipeline.py:1161-1176 plateau_preflight; pipeline.py:1243 env_bad -> refusal() exit 5 (strict); neg test_negative.py:154-156 |
| A2 gamma_GS pre-flight | 1227-1228 | GS/G in [1.4, 1.5] | REFUSE (GateError), every design | **PORTED-GATED** | build | pipeline.py:999-1016 gamma_gs_check (strict refuses, others note) -- only when a meander head is drawn (no GS otherwise); neg test_negative.py:186-199; ALIGNMENT AG7f "partly gated" |
| A3 GS, G >= w_min | 1229-1230 | GS >= w_min and G >= w_min, w_min = max(F_min, h/AR_max) | REFUSE (GateError) | **PORTED-GATED** | build | pipeline.py:1039-1044 wmin_check (G always; GS when a head is drawn); neg test_negative.py:239-248 |
| A4 DIRECT sigma_0 >= 0 | 586-589 | sigma_0 = (1-2f1)/rho1 - f1 >= 0 (w_start not narrower than the closure needs) | REFUSE (GateError) | **EQUIVALENT** | build | floor start: S0 found by bisection upward from 0 (pipeline.py:233-246), so S0 >= 0 by construction; no explicit check |
| A5 lamina overlap | 622-624 | delta*/w_c < 1/2 (S1 path) | REFUSE (WindowError) | **PORTED-GATED** | build | schedule.py:369-374 f_of ValueError -> SpecError (pipeline.py:987-996); trace test_appendix_trace.py:37; neg test_negative.py:168-176 |
| A6 f in (0, 1/2) | 627-628 | 0 < f(i) < 1/2 | REFUSE (WindowError) | **PORTED-GATED** | build | kernels.py:127-128 (x_of_f ValueError); drawn: auditor.py:435-444 P5 `rows_without_cutoff` must be 0; Q2 "f < 1/2" pipeline.py:673-678 |
| A7 M2 row-1 closure tripwire | 633-636 | abs(f(1)(rho(1)+2) - 1) <= eps_id = 1e-10 (FOLD mode) | TRIP (AssertionError) | **EQUIVALENT** | test+build | no in-march tripwire; trace test_appendix_trace.py:139-144 checks 1-P(N) = 2 sigma/(1+2 sigma) < 1e-12 (the same closure identity, at N); auditor P2 mass identity (auditor.py:373-393, hard) |
| A8 M4 PE-closure tripwire | 649-651 | abs(r(w_s)/(rho r(w_c)) - 1) <= 1e3 eps_inv = 1e-9 | TRIP | **EQUIVALENT** | test+build | opencif closes on the exact gap-drop field instead (schedule.py:781-788); trace "network gives every designed f" < 1e-8 (test_appendix_trace.py:127-133); auditor P3 drawn side R vs exact and PE closure <= 0.15 (auditor.py:394-418); Q1 |
| A9 E11 thin lamina | 657-659 | delta(i)/h <= eps_lam = 0.10 on every S1 row | REFUSE (GateError) | **PORTED-GATED** | build | build.py:905-907,932-934 `lamina_ok`, `lamina_all_rows_ok` (cutoff rows and every planned row); pipeline.py:1251-1253; parity with CIF_v2.march test_appendix_trace.py:213-247; neg test_negative.py:217-223,324 |
| A10 M16 quasi-developed taper | 660-661 | abs(w_c(i) - w_c(i-1))/Lambda <= eps_taper = 1e-2, Lambda = G + g_p, w_c(0) = W(0) | REFUSE (GateError) | **PORTED-GATED** | build | build.py:908,935-936 `taper_ok`, `taper_rate_max`; parity test_appendix_trace.py:248-260; neg test_negative.py:224-237,332 |
| A11 march reaches F_target | 672-673 | within max_rows = 400000 | REFUSE (RuntimeError) | **PORTED-GATED** | build | schedule.py:725-726 (400000 rows) -> SpecError (pipeline.py:987-996) |
| A12 M5 mass identity | 676-678 | abs(P + 2S - (1 + 2 sigma_0)) <= 1e-10 | TRIP | **PORTED-GATED** | build+test | auditor P2 (auditor.py:376,385: < 1e-10, role hard, every build); trace every row < 1e-12 (test_appendix_trace.py:93-96) |
| A13 N fixed point converges | 1238-1248 | abs(N_new - N) <= 1 within 200 damped iterations | REFUSE (RuntimeError) | **NOT-APPLICABLE** | - | opencif marches row by row on the row's own w_c (schedule.py:678-726); no W(t = i/N) sampling, hence no fixed point on N |
| A14 M7 inverse round trip in design() | 1252-1256 | max abs(f_B/f_A - 1) <= 1e-6 (sigma-form inverse map) | TRIP | **EQUIVALENT** | build | replaced by an independent network solve of the drawn mask: Q1 (pipeline.py:659-667; every gap within 0.25 %, F within 1e-4); ALIGNMENT AA18 |
| A15 M6 three F forms (sigma_0 = 0) | 1259-1263 | abs(F - 2S(N)) and abs(F - 2/(rho(N)+2)) <= 1e-10 | TRIP | **PORTED-GATED** | test | test_appendix_trace.py:139-144 (< 1e-12) |
| A16 terminal w_s(N) >= w_min | 1264-1265 | w_s(N) >= w_min | REFUSE (GateError) | **PORTED-GATED** | build | pipeline.py:1045-1049 wmin_check (planned); drawn: assess aspect ratio over measured widths incl. outlet leads (pipeline.py:719-744) |
| A17 foot u_min >= u_ff | 1271-1273 | u_min >= 0.012 | REFUSE (WindowError) | **PORTED-GATED** | build | pipeline.py:747 env `foot`, 1249-1250 env_bad (strict refuses); not seeded by a negative test (ALIGNMENT sec. 5) |
| A18 M13 gap dominance | 1279-1283 | max R_gap f/R_c < 0.05 | REFUSE (GateError) | **PORTED-GATED** | build | build.py:862,923 `M13_gap_dominance_max`, `M13_ok`; pipeline.py:1244; seeded only with other gates (ALIGNMENT sec. 5) |
| A19 M14 channel Re | 1284-1288 | max Re <= 10, Re = rho Q_c D_h/(mu h w_c), Q_c = Q_in P(i-1), every row | REFUSE (GateError), always | **PORTED-GATED** | build | build.py:813-814,923 `Re_max`, `Re_ok`; drawn mask: auditor P8-M14 (auditor.py:478,493-504, role envelope); neg test_negative.py:157 |
| A20 M14 Stokes gates | 1289-1300 | Re_gap = Re(f1 Q_in; G x h) <= 1e-2 and Re_delta = rho C(alpha) Q_in/(a b^2) delta^2/mu <= 1e-2, row 1; S1 designs only | REFUSE (GateError) | **PORTED-GATED** | build+test | build.py:818-827,859-861,927 `stokes_ok` (every row, CIF_v2's row-1 literal included); literal parity test_appendix_trace.py:193-211; drawn: auditor P8-M14; neg test_negative.py:158-167 |
| A21 Process: l_g >= g_p (stadium, not bowtie) | 204-210 | l_g >= g_p | REFUSE (GateError, constructor) | **PORTED-GATED** | build | spec.py:599-603 SpecError (post_width_um >= post_length_um); kernels.py:205-208 stadium_post GateError (test_kernels.py:81-87) |
| A22 Pairing.kernel in {cubic, flam} | 242-244 | kernel name | REFUSE (ValueError) | **EQUIVALENT** | build | no kernel option; kernels.py:129-130 x_of_f refuses a non-cubic kernel; trace "pairing rule", "kernel declaration" (test_appendix_trace.py:52-57) |
| A23 r(): positive inputs, bracket > 0 (K2) | 341-346 | w, h, mu > 0; B > 0 | REFUSE (ValueError) | **PORTED-GATED** | build | kernels.py:74-79, bit-identical (test_kernels.py:47); spec.positive_values (spec.py:173-215) |
| A24 r_inverse bracketing | 370-380 | <= 200 halvings/doublings; lo >= 1e-12 | REFUSE (ValueError) | **PORTED-GATED** | build | kernels.py:92-102 (test_kernels.py:51) |
| A25 x_of_f domain | 483-484 | 0 < f < 1/2 | REFUSE (ValueError) | **PORTED-GATED** | build | kernels.py:127-128; test_kernels.py:60-67 |
| A26 route_fold G8 stroke pitch | 870-872 | Lambda >= p = GS + max(F_min, 2 R_min) | REFUSE (RoutabilityError; fold row skipped) | **NOT-APPLICABLE** | - | opencif draws no fold; a row that would fold (straight w_s < w_min) is refused in strict mode (pipeline.py:1053-1063) |
| A27 route_fold G14 own-slot capacity | 896-901 | Phi <= W/p (else arms needed) | REFUSE (RoutabilityError; row skipped, schedule kept) | **NOT-APPLICABLE** | - | as A26 |
| A28 stadium_post bowtie | 2591-2593 | W_perp >= L_along | REFUSE (GateError) | **PORTED-GATED** | build | kernels.py:205-208, used at layout.py:736-749; spec.py:599-603 |
| A29 uturn_widths pull-in | 2639-2641 | 2 (w_s - pull_in) > 0 | REFUSE (GateError) | **EQUIVALENT** | build | layout.loop_profile solves each side for its target R (layout.py:176-216) and fails by name (layout.py:59-65); pull-in is reproduce-only |
| A30 uturn_widths radius | 2642-2645 | inner rail radius R_turn - w_c/2 - l_g - w_s/2 > 0 | REFUSE (GateError) | **EQUIVALENT** | build | layout.py:270 "could not find a loop radius meeting the floor"; build.collisions part_clearance (build.py:981-1014) |
| A31 heleshaw_split validity gate | 3149-3158 | valid iff h < 0.2 min(w) | diagnostic (G19 requires it to fire) | **NOT-APPLICABLE** | - | opencif never computes F by a uniform-sheet law (its 1-D network uses the exact duct law) |
| A32 build_mask port radius cap | 2748-2758 | r_port = min(port_d/2, 0.40 w_channel) so the disc stays in its channel | rendering rule (checked by G25) | **EQUIVALENT** | build | components_drawn: a 0.9 r disc of each outlet port must lie inside the fluid (build.py:227-239); requested_vs_drawn hard (pipeline.py:644-650) |
| A33 transition_seam has a grid | 4890-4891 | pattern has grid units | REFUSE (GateError) | **NOT-APPLICABLE** | - | meander_core grid->stadium seam; opencif has no such seam |
| A34 transition_seam_classic has a grid | 5754-5755 | same | REFUSE (GateError) | **NOT-APPLICABLE** | - | as A33 |
| A35 stadium_post_arcs bowtie | 5537-5538 | W_perp >= L_along | REFUSE (caught: falls back to a polyline) | **NOT-APPLICABLE** | - | true-arc DXF writer; opencif writes polylines only (as CIF_v2's own acceptance emitter export_mask_dxf does) |

## B. Mandatory reports of design() (Device.reports and Device.report())

| ID / name | CIF_v2 line(s) | what it checks, threshold | role in CIF_v2 | opencif 1.0.0 status | level | evidence (file:line, record field) |
|---|---|---|---|---|---|---|
| B1 L_e/Lambda | 1301-1302,1314 | Stokes-floor L_e(Re row 1, D_h(w_c0))/Lambda; "REPORT, do not gate" | REPORT | **PORTED-REPORTED** | build | build.py:863,959 `gates.Le_over_Lambda_L` (kernels.entrance_length bit-identical, test_kernels.py:79) |
| B2 peak wall shear via C(alpha) | 1291-1293,1314 | gamma_max = C(alpha) Q_in/(a b^2) at row 1, "NOT 6" | REPORT | **PORTED-REPORTED** | build | build.py:832-843,953 `gates.shear_peak` (row1_per_s, gamma_max_per_s over all rows, at the design flow) |
| B3 total pressure drop | 1303,1314 | dP = sum R_c(i) Q_c(i) along the trunk, reported in kPa | REPORT | **MISSING** | - | not computed anywhere (no pressure in any record); the ingredient `design.R_element_SI` (gap 1 -> outlets per unit flow, build.py:1311) is recorded |
| B4 gate values | 1313-1314 | gap_dom, Re_max, Re_gap, Re_delta recorded | REPORT | **PORTED-REPORTED** | build | build.py:923-928 (`validation.gates`) |
| B5 D as a band off the plateau | 1274-1276 | u_ff <= u_min < u_lo: D = [kappa_lo delta*, kappa_hi delta*] | REPORT | **EQUIVALENT** | build | build.py:882-904 `gates.cutoff_band_um`: the kappa interval on the plateau; below it App. B C7's Hill band [kappa(u) delta*, kappa_g delta*] (a different band, ALIGNMENT AB5) |
| B6 F_real = eta F | 1310 | eta = 0.903 | REPORT | **MISSING** | - | not applied and not reported (documented: LIMITATIONS 5; Pairing.eta kept, kernels.py:157) |
| B7 transfer bias (S5, P3) | 260-269,1306 | D_true/D_pred = sqrt(mean c2(beta_cal)/c2(beta(i))) over the rows; "REPORTED, NEVER APPLIED" | REPORT (mandatory) | **MISSING** ? | - | not computed; opencif's disposition calls beta_cal "not needed on the cubic pair" (test_appendix_trace.py:344; ALIGNMENT AA12, AG2), which is the glossary's statement about N, while CIF_v2 reports this bias on D for its default cubic pair |
| B8 sensitivity map (S16-e) | 1101-1132,1307 | dlnF/dln w_s(out), w_c(out), h; dlnD/dln w_c, h, f; dlnN/dln kappa = 2; "MANDATORY REPORT" | REPORT (mandatory) | **MISSING** | - | not computed (the record has dD_last/dln R_side,out and dln R_ext only: pipeline.py:811,816) |
| B9 wall porosity phi = G/Lambda | 199-202,1193 | reported because the kernels assume a solid wall | REPORT | **MISSING** | - | not reported (trivially G/(G + g_p) from the spec) |
| B10 fold plan | 824-853,1305 | Phi(i), L, n_turns, arms against the packing cap W/p | REPORT | **NOT-APPLICABLE** | - | no folds; the head's own capacity plan is schedule.py:801-818 (over_cap, arms_needed) |
| B11 DIRECT mass defect 2 sigma_0 | 1181 | reported | REPORT | **PORTED-REPORTED** | build | floor-start note with S0 (pipeline.py:247-248); `design.spec.S0` in every record (build.py:1298) |
| B12 array length vs budget | 1177 | sum Lambda vs n_legs x leg_length | REPORT | **PORTED-REPORTED** | build | build.py:1289-1297 `validation.length_budget`; pipeline.py:881-888 `status.length_budget_G9` (reported, never failed) |
| B13 u range against the plateau | 1175-1176 | min/max u = D/w_c | REPORT | **PORTED-REPORTED** | build | `gates.plateau_u_min` (build.py:960), `plateau_u_min_all_planned_rows` |
| B14 inverse round-trip value | 1314 | rt | REPORT | **NOT-APPLICABLE** | - | inverse map replaced by the network solve (A14); Q1 values recorded instead |

## C. Supplied limits: Tolerances (11), Pairing (9), Process (19), and the auditor's PhysicsSpec

The trace test accounts for every Tolerances/Pairing/Process field (test_appendix_trace.py:330-379; ALIGNMENT sec. 5).  It does not account for the validators, reports and acceptance functions of the other sections.

| class | field | value | CIF_v2 use | opencif disposition | see |
|---|---|---|---|---|---|
| Tolerances | eps_r | 1e-9 | K4, K5 (validate) | not used (solver tolerance; kernels bit-identical) | test_appendix_trace.py:355-357 |
| Tolerances | eps_inv | 1e-12 | M4 tripwire at 1e3 x | not used (exact closure) | idem |
| Tolerances | eps_cub | 1e-3 | none: defined, listed by R2, never read | not used | idem |
| Tolerances | eps_id | 1e-10 | M2, M5, M6 | same value in auditor P2 (auditor.py:385) |  |
| Tolerances | eps_geo | 1e-9 | G1, G2 (folds) | not used |  |
| Tolerances | eps_net | 2e-3 | terminal layer, N2, N4 | not used (Q1 2.5e-3; the auditor uses ladder_physics' 1e-3) |  |
| Tolerances | eps_lam | 0.10 | E11 (march), CU2 | GATED (build.py:47) | A9 |
| Tolerances | eps_taper | 1e-2 | M16 (march) | GATED (build.py:48) | A10 |
| Tolerances | Re_max | 10 | M14 channel | GATED (build.py:35) | A19 |
| Tolerances | Re_stokes | 1e-2 | M14 Stokes | GATED (build.py:36) | A20 |
| Tolerances | gap_dom | 0.05 | M13 | GATED (build.py:34) | A18 |
| Pairing | kernel | cubic | three call sites | consumed; non-cubic refused | A22 |
| Pairing | kappa | 4.9138 | D = kappa delta | consumed (schedule.py:38) |  |
| Pairing | kappa_ci | (4.750, 5.078) | D band | consumed (cutoff band) | B5 |
| Pairing | beta_cal | (0.80, 0.64, 0.53) | transfer bias report; kappa_for | not used | B7 (MISSING report) |
| Pairing | u_lo | 0.035 | plateau guard; D band | GATED | A1 |
| Pairing | u_ff | 0.012 | foot WindowError | GATED | A17 |
| Pairing | eta | 0.903 | F_real report | not applied, not reported (LIMITATIONS 5) | B6 |
| Pairing | f_soft | 5e-4 | D3 WARN | REPORTED | D3 |
| Pairing | gamma_GS | (1.4, 1.5) | design() pre-flight, every design | PARTLY GATED (when a head is drawn) | A2 |
| Process | h | 140 um | depth | spec input |  |
| Process | mu | 1e-3 Pa s | r; Reynolds | fluid section (gates only); schedule at MU | LIMITATIONS 2 |
| Process | rho_f | 1000 | Reynolds | fluid section | LIMITATIONS 2 |
| Process | F_min | 19 um | w_min; fold pitch; G7 | spec input; DRC floor (gated); part clearance |  |
| Process | R_min | 0 | fold pitch p = GS + max(F_min, 2 R_min) | process.min_radius_um: rounding_check (gated) |  |
| Process | AR_max | 10 | w_min = max(F_min, h/AR_max) | GATED (strict) | A3, A16, G7 |
| Process | G | 20.5 um | gap | spec input (w_min gate) | A3 |
| Process | g_p | 25.6 um | post along flow; Lambda = G + g_p | spec input |  |
| Process | l_g | 50.6 um | wall thickness; R_gap; must be >= g_p | spec input; >= g_p gated; not measured on the drawing (LIMITATIONS 12) | A21 |
| Process | c_ent | 1.0 | R_gap entrance allowance | spec input (gap_entrance_c) |  |
| Process | GS | 28 um | fold width; gamma_GS; w_min | meander.channel_um (head only) | A2, A3 |
| Process | W | 1562 um | fold packing budget per side | meander head width; head inside W/2 (gated) | G5 |
| Process | l_stroke | 1562 um | fold stroke | not applicable (no folds) |  |
| Process | L_eq | 0 | bend allowance | superseded: exact bend; c_corner for square turns | ALIGNMENT AG7b |
| Process | n_legs | 1 | length budget | layout.legs |  |
| Process | leg_length | 75.8 mm | length budget (G9) | chip leg length, reported | B12 |
| Process | Q_in | 30 uL/min | Reynolds only | operation.flow_uL_min (gates) |  |
| Process | m_budget | 30 | N7 | not applicable (exact outlet) | N7 |
| Process | L_out | 200 um | run-out (G13) | outlet.stretch_um; exact outlet closure | G13 |
| PhysicsSpec | tol_f_bulk | 5e-3 | P4 | PORTED (auditor.py:67; read from the source by test_auditor.py:52-55) |  |
| PhysicsSpec | tol_delta | 0.02 | P5 | PORTED (auditor.py:68) |  |
| PhysicsSpec | tol_closure | 5e-4 | none: defined, never read (P3 uses a literal 0.15) | not used (P3 ports the literal 0.15) |  |
| PhysicsSpec | eps_net | 1e-3 | P4 terminal window | PORTED (auditor.py:69) |  |
| PhysicsSpec | Re_max | 1.0 | P8; compile_unified Reynolds pre-flight | PORTED as a report (auditor.py:74) | CU4, AUD-P8 |
| PhysicsSpec | Q_in_uLmin | 5.0 | P8 operating point | PORTED (auditor.py:75) |  |
| PhysicsSpec | w_min, wall_min | 6 / 6 um | P9; compiler | process F_min / wall rules | AUD-P9 |
| PhysicsSpec | max_correction_iters | 8 | F7 | not applicable | F7 |

Constants by disposition: consumed 3, gated 8, not used / not applicable 11, partly gated 1, ported to the auditor 5, reported 1, spec input, replaced or reused by another check 18.

## D. The validate() battery: 162 validators (VALIDATORS registry, 1523-2527 and later registrations)

Role "VAL" = a pass condition of `validate()`, which runs on CIF_v2's own fixtures (`_Ctx`: a short FOLD device, its DIRECT twin, `reference_request()`, the embedded demo specs), never on a user's design.  An opencif counterpart in the build is therefore stronger than the original; a counterpart in `tests/` is the like-for-like.

| ID / name | CIF_v2 line(s) | what it checks, threshold | role in CIF_v2 | opencif 1.0.0 status | level | evidence (file:line, record field) |
|---|---|---|---|---|---|---|
| K1 sorting | 1536 | r(w,h) == r(h,w) bit for bit | VAL | **PORTED-GATED** | test | test_kernels.py:47 r = CIF_v2.r bit for bit, both orders of the sides |
| K2 bracket B | 1544 | B > 0, increasing on [1, 60]; B(1) = 0.4217310 (1e-6, 400 terms) | VAL | **PORTED-GATED** | build+test | auditor P1 B(1) +-5e-7 (auditor.py:81,344-352, hard); duct_bracket bit-identical over alpha in [1, 400] (test_kernels.py:41); monotonicity not re-run |
| K3 r decreasing | 1554 | r strictly decreasing over w in [2, 1400] um | VAL | **EQUIVALENT** | test | inherits from bit-identity with CIF_v2.r (test_kernels.py:47); not re-run |
| K4 truncation | 1563 | abs(r(50 odd)/r(100 odd) - 1) < eps_r = 1e-9 | VAL | **EQUIVALENT** | test | trace: 50-term r against a 401-term series < 1e-6 (test_appendix_trace.py:72; ALIGNMENT AA26: 1.9e-9) |
| K5 symmetry | 1575 | r(w;h) = r(h;w) < 1e-9 | VAL | **PORTED-GATED** | test | as K1 |
| K6 square duct | 1583 | r h^4/mu = 28.454154 (1e-4) | VAL | **PORTED-GATED** | build+test | auditor P1 square_duct_R 28.454154 +-5e-6 (auditor.py:81,345); trace +-5e-7 (test_appendix_trace.py:74) |
| K7 slender asymptote | 1592 | alpha (1 - B(alpha)) -> 0.6302489 at alpha 1e4 (400 terms) | VAL | **EQUIVALENT** | - | same code as CIF_v2.B (test_kernels.py:41); the alpha = 1e4 limit is not re-run |
| K8 inversion round trip | 1600 | abs(r(r^-1(y))/y - 1) < 1e-9 | VAL | **EQUIVALENT** | test | r_inverse bit-identical (test_kernels.py:51); not re-run |
| K9 reciprocity | 1610 | B(1/a) = a^-2 B(a) < 1e-7 | VAL | **EQUIVALENT** | - | code identity; not re-run |
| K10 square-duct f Re | 1621 | 24/B(1) = 56.908 (0.01) | VAL | **PORTED-GATED** | build | auditor P1 12/B(1) = 28.454154 is the same anchor (f Re/2) |
| K11 plate f Re | 1629 | f Re -> 96 (0.3) at alpha 4000 | VAL | **EQUIVALENT** | - | code identity; not re-run |
| E1 F_lam -> cubic | 1639 | F_lam(x, beta -> inf) = 3x^2 - 2x^3 (1e-5) | VAL | **NOT-APPLICABLE** | - | opencif uses the cubic kernel only (ALIGNMENT AA12) |
| E2 hard bound | 1648 | F_lam(1/2, beta) = 1/2 (1e-14) | VAL | **EQUIVALENT** | build | cubic analogue g(1/2) = 1/2 in auditor P1 (auditor.py:82,348) |
| E3 F_lam monotone, A_n > 0 | 1656 |  | VAL | **NOT-APPLICABLE** | - | F_lam not used |
| E4 c2 limits | 1668 | c2(inf) = 3, c2(1) = 3.5568, decreasing | VAL | **NOT-APPLICABLE** | - | c2 not used |
| E5 c2 analytic tail | 1679 |  | VAL | **NOT-APPLICABLE** | - | c2 not used |
| E6 cancellation guard | 1692 |  | VAL | **NOT-APPLICABLE** | - | F_lam not used |
| E7 K_depth | 1708 | 0 < K_depth <= 1 | VAL | **NOT-APPLICABLE** | - | not used |
| E8 inversion seed | 1718 | sqrt(f/c2) within 5 % | VAL | **NOT-APPLICABLE** | - | not used |
| E9 one duct series | 1731 | B is the single call site for r, C | VAL | **EQUIVALENT** | test | kernels.C_shear calls duct_bracket (kernels.py:170-180); bit-identical (test_kernels.py:77) |
| E10 shallow scaling | 1740 | beta c2 -> ~1.7 | VAL | **NOT-APPLICABLE** | - | not used |
| E11 thin lamina on the fixture | 1748 | max delta/h <= 0.10 | VAL | **PORTED-GATED** | build | A9 |
| C1 g monotone | 1759 | g'(x) = 6x(1-x) >= 0 | VAL | **EQUIVALENT** | test | x_of_f bit-identical, a bracketed root on [0, 1/2] (kernels.py:120-131; test_kernels.py:58) |
| C2 g(1/2) = 1/2 | 1767 | to machine precision | VAL | **PORTED-GATED** | build | auditor P1 g_half (auditor.py:82,348, hard) |
| C3 unique root | 1774 | x(f) in (0, 1/2), abs(g(x) - f) < 1e-12 for f in [1e-5, 0.49] | VAL | **PORTED-GATED** | build+test | auditor P1 round trip f in {1e-6 ... 0.49} < 1e-12 (auditor.py:346); trace test_appendix_trace.py:40,44 |
| C4 seed error law | 1783 | abs(x0/x - 1) = sqrt(f/27) within 5 % | VAL | **EQUIVALENT** | test | trace: leading order sqrt(f/3) (1e-3 at f 1e-6) and <= 0.47 % to the cap (test_appendix_trace.py:41,47); opencif never uses a seed |
| P1 one Pairing, three call sites | 1797 | f_of_x / x_of_f round trip; D = kappa delta | VAL | **PORTED-GATED** | test+build | trace test_appendix_trace.py:52-57; `design.kernel_pair` in every record (build.py:66-69,1298) |
| P2 kappa tagged with its kernel | 1808 | kappa_flam = kappa/K_depth(beta_cal) > kappa | VAL | **NOT-APPLICABLE** | - | kappa_F never used (ALIGNMENT AG2) |
| P3 transfer bias reported, never applied | 1818 | band computed and > 0; D uncorrected | VAL | **MISSING** ? | - | D is uncorrected (true) but the bias is not reported: B7 |
| M1 sigma_0 = 0 in FOLD mode | 1830 | bit-exact 0.0 | VAL | **EQUIVALENT** | build | sealed start: Spec.S0 = 0.0 (schedule.py:535), mass identity with S0 (auditor.py:374-376) |
| M2 row-1 closure | 1837 | abs(f(1)(rho(1)+2) - 1) < 1e-10 | VAL | **EQUIVALENT** | test | A7 |
| M3 f in (0,1/2); sigma increasing | 1846 |  | VAL | **PORTED-GATED** | build | A6 (sigma increases whenever f > 0) |
| M4 closure identity on every row | 1855 | < 1e-9 | VAL | **EQUIVALENT** | test+build | A8 |
| M5 mass identity | 1865 | < 1e-10 | VAL | **PORTED-GATED** | build+test | A12 |
| M6 F three forms | 1874 | < 1e-10 | VAL | **PORTED-GATED** | test | A15 |
| M7 round trip, both modes | 1885 | < 1e-9 | VAL | **EQUIVALENT** | build | A14 |
| M8 delta(i) = delta* | 1900 | both kernels, < 1e-9 | VAL | **PORTED-GATED** | test+build | cubic: test_appendix_trace.py:100-102; drawn: Q2 cutoff at every gap (pipeline.py:671-678) |
| M9 invariance sweeps | 1914 | F, N unchanged under mu x2, Q_in x2, g_p (Lambda), GS | VAL | **EQUIVALENT** ? | - | by construction: schedule.build never reads Q_in, mu is the constant MU (schedule.py:36), f is assigned from delta*/w_c; no invariance test |
| M10 fixed point, non-monotone W | 1944 |  | VAL | **NOT-APPLICABLE** | - | no arbitrary W(t) (tapers: constant, constant shear; LIMITATIONS 16) |
| M11 N a true fixed point | 1957 | seed x3 -> abs(dN) <= 1 | VAL | **NOT-APPLICABLE** | - | no fixed point (A13) |
| M12 Psi canonical taper | 1969 | Psi linear in i; F_Psi vs march < 1e-3 | VAL | **NOT-APPLICABLE** | - | Psi taper not offered (LIMITATIONS 16; ALIGNMENT AA25) |
| M13 gap dominance (and it fires) | 2000 | < 0.05 | VAL | **PORTED-GATED** | build | A18 |
| M14 Stokes gates on the fixture | 2016 | Re_gap, Re_delta < 1e-2; Re < 10 | VAL | **PORTED-GATED** | build | A19, A20 |
| M15 entrance-length Stokes floor | 2026 | L_e/D_h -> 0.619 as Re -> 0 | VAL | **PORTED-REPORTED** | build+test | B1; bit-identical entrance_length (test_kernels.py:79) |
| M16 taper rate on the fixture | 2038 | <= 1e-2 | VAL | **PORTED-GATED** | build | A10 |
| M17 sensitivity map emitted | 2048 | keys present; abs(dlnF/dln h) < abs(dlnF/dln w_s)/3 with both outlet channels < h | VAL | **MISSING** | - | B8 |
| M18 transport lemma | 2063 | abs(2 phi_crit - F) < 1e-6 | VAL | **PORTED-GATED** | test | test_appendix_trace.py:97-99 (CIF_v2.critical_phi on opencif's schedule) |
| M19 Gifford Fig. 2, engine-driven | 3368 | N vs i_T <= 1 row; Q_s/Q_c vs [(1-2f)^-i - 1]/2 < 5e-3; depth collapse < 5e-3; w_s fans > 5 % | VAL | **EQUIVALENT** ? | test | partial: the implementers' constant-f march (f 1e-4: P, sigma, F at row 5000; test_appendix_trace.py:76-87) checks the closed form at one row; no i_T or depth-collapse check |
| D1 u >= u_ff, and it fires | 2077 |  | VAL | **PORTED-GATED** | build | A17 |
| D2 band off the plateau | 2090 |  | VAL | **EQUIVALENT** | build | B5 |
| D3 f <= f_soft, else WARN | 2101 | f_soft = 5e-4 | VAL | **PORTED-REPORTED** | build | build.py:931 `f_above_soft`; note pipeline.py:773-775 |
| D4 delta <= w_c/2 | 2111 |  | VAL | **PORTED-GATED** | build | A5 |
| G1 resistance DRC on the drawn traversal | 2120 | abs(r(GS) L_measured/(rho r(w_c) Lambda) - 1) < 1e-6 (routed folds) | VAL | **EQUIVALENT** | build | the folded side channels opencif draws are head units: each unit's length and width measured on the polygon (measure.py:91ff), `fidelity.unit_length_rel_max` reported (build.py:1231-1239), Q1 gated |
| G2 series topology (no shortcut) | 2143 | L >= Phi Lambda (1 - eps_geo) and L > 5 x straight distance | VAL | **EQUIVALENT** | build | topology: 2 holes per head unit expected and found (build.py:256-263; pipeline.py:625-626); features_found; Q1 |
| G3 no stub | 2156 | fold starts at link i, ends at link i+1 | VAL | **EQUIVALENT** | build | as G2 |
| G4 constant cross-section | 2167 | axis-aligned strokes => width GS everywhere | VAL | **EQUIVALENT** | build | unit widths measured by cross-sections; where G4 fails (square corners) the corner term is computed (LIMITATIONS 8) |
| G5 packing | 2178 | arms == 0 iff Phi <= W/p; arms counted | VAL | **EQUIVALENT** | build | head capacity: schedule.py:801-818 (over_cap, arms_needed, A_max); head inside W/2 (pipeline.py:651, hard) |
| G6 global area (necessary only) | 2191 | sum Phi Lambda vs (W/p) sum Lambda; informational | VAL (informational) | **NOT-APPLICABLE** | - | fold area bound |
| G7 minimum feature | 2201 | every w_c >= w_min; w_s(N) >= w_min; GS, G >= w_min; g_p >= F_min | VAL | **PORTED-GATED** | build | pipeline.py:1028-1063 wmin_check (G, GS, w_c, w_s(N), every straight w_s); drawn: aspect ratio (pipeline.py:719-744), DRC solid at F_min - 0.1 (posts, walls) and auditor P9 (auditor.py:516-528) |
| G8 fold pitch fires | 2212 | Lambda >= p | VAL | **NOT-APPLICABLE** | - | A26 |
| G9 length budget | 2228 | array length vs n_legs x leg_length; reports OK/OVER | VAL (report only) | **PORTED-REPORTED** | build | B12 |
| G10 turn rebalance with the upstream rho | 2239 | r(w_in) L_in = r(w_out) L_out = rho r(w_c) L_p (< 1e-9), w_in < w_out | VAL | **EQUIVALENT** | build | layout.loop_profile: each side's integral of r(w) along its own centreline = the row's exact target R_s(k) (layout.py:176-216), with the exact curved-duct correction |
| G11 handover: folded iff w_s < w_min | 2254 |  | VAL | **EQUIVALENT** | build | predicate kept, consequence changed: a straight w_s < w_min is refused (strict) instead of folded (pipeline.py:1053-1063; LIMITATIONS 7) |
| G12 bend accounting | 2265 | n_turns >= 0 from the drawn schedule; L_drawn <= L_ideal | VAL | **EQUIVALENT** | build | corners recorded per unit (geometry.unit_path; LIMITATIONS 8: G12 holds by construction); exact bend solution |
| G13 run-out match | 2273 | (R_out_s/R_out_c)/rho(N) = 1 (1e-9); 3x mismatch moves F > 30 % | VAL | **EQUIVALENT** | build | exact outlet closure sizes the leads (layout.py:1280-1310); auditor P7 x3 sensitivity (report); model scenario outlet +-1 % (pipeline.py:802-811) |
| G14 unroutable fold raises | 2295 |  | VAL | **NOT-APPLICABLE** | - | A27 |
| N1 Laplacian solves | 2311 | finite, 0 < F_em < 1 | VAL | **PORTED-GATED** | build | measure.network on every build (build.py:1228); auditor KCL residual < 1e-9 (auditor.py:385) |
| N2 bulk agreement | 2320 | max bulk abs(f_em/f - 1) < eps_net = 2e-3 | VAL | **PORTED-GATED** | build | Q1 at every gap (2.5e-3 on the cutoff; pipeline.py:661); auditor P4 5e-3 vs the kernel |
| N3 F_em = F | 2338 | rel < max(2e-2, 2 m/N) | VAL | **PORTED-GATED** | build | Q1 F within 1e-4 of intended; Q2 F within F_tol (default 1e-3); auditor P6 |
| N4 M13 <-> network coupling | 2352 | gap_dom < 0.05 => bulk < eps_net; violated => bulk moves | VAL | **NOT-APPLICABLE** | - | self-test of CIF_v2's PE closure; opencif's exact closure agrees with its network independently of M13 (trace < 1e-8) and Q1 checks every mask; the mouth R_gap +-20 % scenario is reported (scenarios.py:64-69) |
| N5 Q_s(0) = 0 structural | 2376 | f_em[0] > 0 with no upstream side edge | VAL | **EQUIVALENT** | build | sealed start: the side rails begin at gap 1 in measure.solve3 and auditor.solve_ladder (auditor.py:164-189) |
| N7 terminal layer within budget | 2385 | m_meas <= m_budget = 30 | VAL | **NOT-APPLICABLE** | - | exact outlet closure (disposition of m_budget, test_appendix_trace.py:369); auditor n_term reported (auditor.py:431-432) |
| N8 symmetry precondition | 2395 | G10 registered | VAL | **EQUIVALENT** | build | both sides solved explicitly, no symmetry assumed (measure.py docstring; auditor.solve_ladder) |
| N9 boundary amplitude; q non-monotone | 2402 | report | VAL (report) | **PORTED-REPORTED** | build | auditor terminal_window q_min, lambda, n_term in P4 detail (auditor.py:192-199,431-432) |
| S1 C(alpha) -> 6 | 2420 | C(1e6) = 6 (1e-3) | VAL | **EQUIVALENT** | test | C_shear bit-identical (test_kernels.py:77); not re-run |
| S2 C(1), C(2) | 2428 | 9.6078, 8.1341; 6 understates by 26.2 % | VAL | **EQUIVALENT** | test | as S1 |
| S3 C shares B | 2438 |  | VAL | **EQUIVALENT** | test | kernels.py:170-180 |
| S4 constant Q_gap is not constant shear | 2449 | C drifts > 5 % over w_c 75 -> 52 | VAL | **PORTED-REPORTED** | build | `shear_peak.over_the_taper` max_over_min (build.py:844-858; ALIGNMENT S1) |
| R1 determinism | 2462 | same inputs -> same N, F | VAL | **EQUIVALENT** ? | - | rebuild comparisons in provenance (byte-identical photomasks, CHANGELOG rc8); no determinism test in tests/ |
| R2 tolerances explicit | 2471 | Tolerances has the named fields | VAL | **PORTED-GATED** | test | test_kernels.py:89-91; trace dispositions test_appendix_trace.py:350-357 |
| R3 constants tagged and echoed | 2480 | kernel, kappa, beta_cal echoed; UNCHECKABLE >= 8 | VAL | **EQUIVALENT** | build+test | `design.kernel_pair` in every record; PRESET_PROVENANCE codes external inputs A; trace field accounting |
| G15 stadium posts | 3172 | simple polygon, area between ellipse and rectangle; 2 posts per row | VAL | **EQUIVALENT** | build+test | stadium_post bit-identical (test_kernels.py:81); 2 holes per straight row expected and found (build.py:256-263, hard) |
| G16 equal-resistance U-turn rails | 3188 | R_in/R_out = 1 (1e-9); linear-in-radius rule visibly unbalanced | VAL | **EQUIVALENT** | build | G10 |
| G17 connectivity | 3212 | largest component > 0.999 of open pixels | VAL | **PORTED-GATED** | build | DRC n_components == 1 and topology single polygon (pipeline.py:617-626, hard) |
| G18 geometric F-test | 3226 | widths measured off the polygons: dwc, dws < 2e-2; abs(F_poly - F_design) < 5e-2 | VAL | **PORTED-GATED** | build | measure_device + network on the drawn polygon, Q1 (every gap 2.5e-3, F 1e-4; pipeline.py:659-667) |
| G19 Hele-Shaw gate fires | 3242 | invalid, error > 5 % | VAL | **NOT-APPLICABLE** | - | A31 |
| G20 packing bound; min_w_start lever | 3256 | Phi(1) at min_w_start = W/p (1e-9), Lambda-free | VAL | **NOT-APPLICABLE** | - | fold packing |
| G21 composition seams | 3275 | one march, 3 renderings: d_D < 1e-9, seam step/neighbour median < 3 | VAL | **NOT-APPLICABLE** | - | series stages are not offered (LIMITATIONS 17); the one seam opencif draws (head -> linear) is inside one march and covered by Q2 at every gap |
| G22 DXF round trip | 3389 | LWPOLYLINE and vertex counts, AC1015, $INSUNITS | VAL | **PORTED-GATED** | build | auditor P0: the file read back equals the drawn polygon vertex for vertex, units 13, closed polylines, no arcs (auditor.py:333-342, hard) |
| G23 architecture matrix | 3411 | linear, serpentine 4/8, paperclip 5 and composed: build, connect > 0.999, F from polygons < 5e-2, DXF, 4 ports | VAL | **EQUIVALENT** | build+test | every build runs topology, DRC, Q1; the layouts are exercised by test_rc7.py, test_composition.py, the presets |
| G24 depth-correct field solve | 3606 | F_field within 25 % of design; uniform sigma > 3x worse | VAL | **EQUIVALENT** ? | build | no 2-D field solve on the mask; a short across a wall is caught by the hole count (a cut post changes it), the per-gap width measurement and Q1 |
| G25 ports inside their channels | 3633 | INLET, OUT_C, OUT_S+, OUT_S-; each disc > 85 % open | VAL | **EQUIVALENT** | build | A32; port_separation (build.py:266-287, hard); combined outlets have one filtrate port by design |
| G26 wall integrity | 3656 | wall-band open fraction = link + 4 cap corners within 5 % | VAL | **EQUIVALENT** | build | hole count, DRC solid band at F_min, features_found (pipeline.py:653-658), gap widths measured at mid-wall |
| X1 embedded kernels = CIF kernels | 3770 | B 1e-9, F_lam 1e-5, c2 5e-6, C 5e-5 | VAL | **EQUIVALENT** | test | test_kernels.py: opencif kernels bit-identical to CIF_v2's (B, r, r_inverse, g, x_of_f, C_shear, entrance_length) |
| X2 two marches agree (F_lam) | 3785 | abs(dN) <= 2, dF < 2e-3; doc-form error flagged | VAL | **NOT-APPLICABLE** | - | F_lam / Psi not used |
| X3 v5 anchors in the embedded core | 3810 | B(1), c2(inf), F_lam(1/2), C(1), f Re | VAL | **EQUIVALENT** | build+test | B(1), square duct, g(1/2) in auditor P1; C bit-identical |
| Z1 P1-P10 on demo A all pass | 3835 | every P check ok | VAL | **EQUIVALENT** | test+build | test_auditor.py:123-170 runs CIF_v2's audit() on demos A, B, C and requires the port's numbers (network, window, P3 PE part, P4-P8); opencif's auditor runs on every photomask (build.py:1334) |
| Z2 P1-P10 on demo B | 3854 | every P check ok | VAL | **EQUIVALENT** | test | as Z1 |
| Z3 embedded run_selftest passes | 3869 |  | VAL | **EQUIVALENT** | test | test_auditor.py:233-303 seeded faults (P0, P1, P3, P4/P5, P6, P8-M14, P9, P9-AR, P10) |
| Y1 meander check_params_feasible | 3891-3910 (Y), _MEANDER_SRC 1211 | derive() accepts the Params and reproduces 13 derived fields (1e-9) | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | meander_core (traced line-drawing pattern) geometry; opencif's head is its own (shapes.py, layout.py) |
| Y2 meander check_counts_and_kinds | 3891-3910 (Y), _MEANDER_SRC 1229 | A arms, T taper ... units per side, densely indexed, 7/5/3 points by kind | VAL (on demo A's compiled pattern) | **EQUIVALENT** | build | topology row_kinds_survive (plan kinds == drawn kinds) and holes_expected (build.py:1183-1188) |
| Y3 meander check_manhattan | 3891-3910 (Y), _MEANDER_SRC 1282 | every segment axis-aligned, non-degenerate | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | opencif heads have arcs or square corners with roundings |
| Y4 meander check_row_grid | 3891-3910 (Y), _MEANDER_SRC 1298 | rows on the declared spacing; slot - width >= min_clear | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | meander_core (traced line-drawing pattern) geometry; opencif's head is its own (shapes.py, layout.py) |
| Y5 meander check_junctions_on_trunk | 3891-3910 (Y), _MEANDER_SRC 1358 | every out line starts on the trunk face (1e-6) | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | meander_core (traced line-drawing pattern) geometry; opencif's head is its own (shapes.py, layout.py) |
| Y6 meander check_arm_shape | 3891-3910 (Y), _MEANDER_SRC 1379 | nested arms: reach step = fan_step, white >= min_clear, last rise >= 1.5 pitch | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | meander_core (traced line-drawing pattern) geometry; opencif's head is its own (shapes.py, layout.py) |
| Y7 meander check_taper_shape | 3891-3910 (Y), _MEANDER_SRC 1438 | taper reaches from the profile, strictly decreasing | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | meander_core (traced line-drawing pattern) geometry; opencif's head is its own (shapes.py, layout.py) |
| Y8 meander check_decay_shape | 3891-3910 (Y), _MEANDER_SRC 1490 | exponential decay law re-derived | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | meander_core (traced line-drawing pattern) geometry; opencif's head is its own (shapes.py, layout.py) |
| Y9 meander check_release_shape | 3891-3910 (Y), _MEANDER_SRC 1542 | release units and bridges | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | meander_core (traced line-drawing pattern) geometry; opencif's head is its own (shapes.py, layout.py) |
| Y10 meander check_morph_shape | 3891-3910 (Y), _MEANDER_SRC 1611 | morph schedule non-increasing, fusion band law | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | meander_core (traced line-drawing pattern) geometry; opencif's head is its own (shapes.py, layout.py) |
| Y11 meander check_grid_shape | 3891-3910 (Y), _MEANDER_SRC 1692 | grid rung + rail units | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | meander_core (traced line-drawing pattern) geometry; opencif's head is its own (shapes.py, layout.py) |
| Y12 meander check_chain_links | 3891-3910 (Y), _MEANDER_SRC 1744 | stop-short = merge_off, white pocket >= 1 um, jog lands interior | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | meander_core (traced line-drawing pattern) geometry; opencif's head is its own (shapes.py, layout.py) |
| Y13 meander check_port_interface | 3891-3910 (Y), _MEANDER_SRC 1820 | port cross-section; boundary_w >= min_clear; rails >= 2 pitch | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | meander_core (traced line-drawing pattern) geometry; opencif's head is its own (shapes.py, layout.py) |
| Y14 meander check_symmetry | 3891-3910 (Y), _MEANDER_SRC 1881 | left = mirror(right) (1e-6) | VAL (on demo A's compiled pattern) | **NOT-APPLICABLE** | - | opencif measures and solves both sides separately, assuming no symmetry |
| Y15 meander check_clearance | 3891-3910 (Y), _MEANDER_SRC 1901 | foreign strokes keep >= min_clear (= line_w) of white | VAL (on demo A's compiled pattern) | **EQUIVALENT** | build | DRC solid band at F_min - 0.1 (drc.py:123-127) and part_clearance >= F_min (layout.py:1576) |
| Y16 meander check_inside_box | 3891-3910 (Y), _MEANDER_SRC 1983 | ink inside the margin box | VAL (on demo A's compiled pattern) | **EQUIVALENT** | build | head inside W/2 (pipeline.py:651, hard); chip area reported (pipeline.py:702-711) |
| Y17 meander check_connectivity | 3891-3910 (Y), _MEANDER_SRC 2014 | one component; contact count = 2 (2U - 1 - R) | VAL (on demo A's compiled pattern) | **EQUIVALENT** | build | topology single polygon with the expected holes; DRC one component |
| Y18 meander check_determinism | 3891-3910 (Y), _MEANDER_SRC 2086 | rebuild reproduces geometry and SVG bytes | VAL (on demo A's compiled pattern) | **EQUIVALENT** | build | provenance rebuild comparisons (no test); records carry the photomask sha256 (auditor.py:549) |
| Y19 meander check_width_profile | 3891-3910 (Y), _MEANDER_SRC 2106 | stroke widths = declared profile (1e-6); trunk >= 2 max(min_clear, line_w/2) | VAL (on demo A's compiled pattern) | **EQUIVALENT** | build | widths measured on the drawing (fidelity, build.py:1239-1245; Q1) |
| U1 every meander mutation caught | 3921 | 11 MUTATIONS each flip >= 1 meander check | VAL | **EQUIVALENT** | test | opencif's own teeth: test_auditor.py faults, test_negative.py seeded defects, test_units.py:101-107 |
| U2 targeted mutation groups | 3945 | morph 2, grid 2, release 2 | VAL | **NOT-APPLICABLE** | - | meander_core regions |
| U3 corner radius cap | 3970 | corner_r above 0.45 x smallest white feature refused; fillets rendered with round joins/caps | VAL | **EQUIVALENT** | build | spec.py:545-567 refuses an explicit fillet at or above half the narrowest drawn feature; rounding_check (build.py:119-123, hard) |
| W1 function-modifiable profiles | 4009 |  | VAL | **NOT-APPLICABLE** | - | no user profiles (LIMITATIONS 16) |
| W2 morph transitions | 4037 |  | VAL | **NOT-APPLICABLE** | - | meander_core |
| W3 compactification | 4061 | decay units lever N; grid-dominant demo passes audit | VAL | **NOT-APPLICABLE** | - | meander_core |
| W4 plain CIF two ways | 4088 | straight-section mask connects > 0.999; grid-dominant demo passes audit | VAL | **EQUIVALENT** | build | head-less designs (e.g. examples/gifford_physics_linear.yaml) pass topology, DRC and the audit |
| T1 unified compiler | 4269 | gates pass, P checks pass, schedule = physics_core march | VAL | **NOT-APPLICABLE** | - | CIF_v2's two-compiler reconciliation |
| T2 plateau gate refuses demo C | 4287 | w_c0 75, D 2.17: WindowError; builds with gates off | VAL | **PORTED-GATED** | build+test | the same lineage refusal: plateau gate (A1); Dinh/Xia refused in strict (ALIGNMENT AB20) |
| T3 UnifiedDevice carries both worlds | 4310 | pattern passes the 19 checks; boundary_w = b_wall | VAL | **NOT-APPLICABLE** | - | meander_core |
| V1 corner_r leaves N, F invariant | 4448 |  | VAL | **NOT-APPLICABLE** | - | MeanderSpec API |
| V2 standalone Params surface | 4472 |  | VAL | **NOT-APPLICABLE** | - | MeanderSpec API |
| V3 physics-determined fields protected | 4498 | MeanderSpec exposes no side/centre/taper profile | VAL | **EQUIVALENT** | build | spec.SCHEMA: f_gap is refused (test_appendix_trace.py:152-157); reproduction-only fields need gates.mode reproduce (spec.py:82,590-591) |
| R2a P10 on a compiled device | 4680 | corner_r within cap; fillets; port squared | VAL | **NOT-APPLICABLE** | - | meander_core rendering (opencif's P10 is a different test: AUD-P10) |
| R2b grid has no stadium posts | 4699 |  | VAL | **NOT-APPLICABLE** | - | meander_core |
| R2c rounded raster connected and different | 4721 | > 0.999; > 50 px differ | VAL | **NOT-APPLICABLE** | - | meander_core raster |
| R3a fillet discs tangent | 4787 | < 1e-9 | VAL | **NOT-APPLICABLE** | - | meander_core |
| R3b anti-aliasing converges | 4812 |  | VAL | **NOT-APPLICABLE** | - | raster |
| R3c ink area converges | 4839 |  | VAL | **NOT-APPLICABLE** | - | raster |
| R4a binary grid->stadium seam | 5051 |  | VAL | **NOT-APPLICABLE** | - | seam |
| R4b posts are islands | 5069 | intrusion <= 2 %; carve < 60 % | VAL | **NOT-APPLICABLE** | - | seam |
| R4c fluid crosses through the gaps | 5122 | > 0.95 connected below the seam | VAL | **NOT-APPLICABLE** | - | seam |
| R4d stadium orientation | 5154 | across/along >= 1 | VAL | **NOT-APPLICABLE** | - | seam (opencif: A21, A28) |
| Q1 (CIF_v2) demo hits its targets | 5185 | abs(F_em - F_target) < 0.02; bulk abs(delta_em/delta* - 1) < 5e-3; P5, P6 pass | VAL | **PORTED-GATED** | build | opencif Q2: F within F_tol (1e-3) and every undeclared gap within +-1 % (pipeline.py:668-678); auditor P5, P6 hard |
| Q2 (CIF_v2) fabrication bridge physics-neutral | 5212 | stroke widths > 0; trunk = compiled (1e-6) | VAL | **EQUIVALENT** | build | opencif Q1 (drawn vs intended network at every gap) |
| Q3 (CIF_v2) seam preserves link resistance | 5236 | R_gap/R_rung within 2 % | VAL | **NOT-APPLICABLE** | - | seam |
| R5a meander DXF true arcs | 5445 | CIRCLE count = joins + fillets | VAL | **NOT-APPLICABLE** | - | opencif writes straight-segment polylines only; P0 requires no arcs (auditor.py:335) -- as CIF_v2's acceptance emitter export_mask_dxf (7356) does |
| R5b DXF faithful to the raster | 5470 | connected; ink fraction within 0.8 % | VAL | **EQUIVALENT** | build | P0 exact equality (stronger) |
| R5c DXF valid R2000, um units | 5495 |  | VAL | **PORTED-GATED** | build | P0 units 13, closed polylines, valid polygon (the R2000 version string is not checked) |
| R6a classic DXF stadium caps as ARCs | 5652 | 2 arcs + 2 lines per post | VAL | **NOT-APPLICABLE** | - | as R5a |
| R6b DXF posts = the 2:1 pill | 5676 | g_p and l_g read back within 0.5 um; ratio within 0.05 | VAL | **EQUIVALENT** ? | build | by construction (posts drawn by stadium_post(g_p, l_g), layout.py:736-749; P0 file = drawing; G measured at mid-wall); l_g itself is not measured (LIMITATIONS 12) |
| R6c classic DXF valid | 5703 |  | VAL | **PORTED-GATED** | build | P0 |
| R7a classic-wall 2:1 posts | 5838 |  | VAL | **NOT-APPLICABLE** | - | seam |
| R7b wall step preserves R_gap | 5859 | within 2 % | VAL | **NOT-APPLICABLE** | - | seam |
| R7c classic section connected | 5878 | > 0.95 | VAL | **NOT-APPLICABLE** | - | seam |
| R8a Gifford device does not taper | 5920 | constant w_c | VAL | **EQUIVALENT** | build | constant-taper presets (centre_channel.taper constant) |
| R8b straight side channels, posts in the walls | 5934 | 2 continuous side channels; connected > 0.99 | VAL | **PORTED-GATED** | build | side_continuity: one strip per leg and side, narrowest >= F_min (build.py:290-355; pipeline.py:627-632, hard) |
| R8c serpentine + paperclip join via U-turns | 5961 | connected > 0.95; arc DXF | VAL | **EQUIVALENT** | build | topology single polygon on every layout; arcs N/A |
| Q4 (CIF_v2) F from drawn Gifford polygons | 6037 | dwc, dws < 1e-3; abs(F_poly - F_design) < 0.02 | VAL | **PORTED-GATED** | build | opencif Q1 |
| Q5 (CIF_v2) F from the drawn meander | 6056 | abs(F_em - F_target) < 0.02; P4 pass | VAL | **EQUIVALENT** | build | auditor P4, P6 on every photomask |
| Q6 (CIF_v2) seam gap measured from pixels | 6076 | R_gap/R_rung within 15 % | VAL | **NOT-APPLICABLE** | - | seam |
| R9a serpentine U-turns alternate | 6194 | right, left, right; connected > 0.99 | VAL | **EQUIVALENT** | build | layout construction; a crossing turn fails part_clearance (refusal); topology |
| R9b inlet + three outlets | 6220 | INLET, OUT_C, two OUT_S | VAL | **EQUIVALENT** | build | components_drawn outlet_ports; requested_vs_drawn; port_separation (pipeline.py:633-650) |
| R9c silhouette IoU vs reference PNG | 6240 | IoU > 0.30 (skipped if the PNG is absent) | VAL | **EQUIVALENT** | test | test_trace_anchors.py, test_figure_anchors.py against the papers' figures (stronger) |

## E. The v2 integration layer: unified-compiler gates, build_mask_v2, J-suite, D-suite, optimizer criterion, G14 routability, the V_* acceptance validators, run_acceptance

run_acceptance's pass condition is `ok = goldens and G14 and jsuite` (8041-8043); the D-suite, V_operating_regime, V_steric and V_fab_tolerance are stored in the acceptance JSON but are not part of `ok`.

| ID / name | CIF_v2 line(s) | what it checks, threshold | role in CIF_v2 | opencif 1.0.0 status | level | evidence (file:line, record field) |
|---|---|---|---|---|---|---|
| CU1 compile_unified plateau | 4178-4189 | max w_c <= D/u_lo | PRE (WindowError) | **PORTED-GATED** | build | A1 |
| CU2 compile_unified thin lamina | 4191-4198 | delta*/h <= eps_lam (+1e-12) | PRE (GateError) | **PORTED-GATED** | build | A9 |
| CU3 compile_unified closure == network | 4200-4215 | audit P3 and P4 must pass | PRE (GateError) | **PORTED-GATED** | build | auditor P3, P4 hard: a FAIL (exit 2) after drawing rather than a refusal before it (auditor.py:394-432; pipeline.py:686-693) |
| CU4 compile_unified Reynolds | 4217-4230 | inlet Re <= spec.Re_max = 1.0 at Q_in_uLmin = 5 (MeanderSpec defaults 4397-4400) | PRE (GateError) | **PORTED-REPORTED** | build | deliberately reported, not gated: gates.P8 (build.py:911-918,937) and auditor P8 (auditor.py:505-515); LIMITATIONS 6 |
| MV1 build_mask_v2 no-floor guard | 6571-6580 | side_mode straight: every w_s >= GS unless allow_floor | REFUSE (GateError) | **EQUIVALENT** | build | opencif never floors a drawn width: it draws the planned width and refuses w_s < w_min in strict mode (pipeline.py:1053-1063); drawn widths measured (fidelity, Q1) |
| MV2 composition interference gate | 6483-6543,6707-6749 | turn annuli, turn vs non-adjacent leg band, leads vs turns: clearance >= F_min; mirror/widening search up to +5 mm, else GateError | REFUSE (GateError) | **PORTED-GATED** | build | build.collisions/part_clearance: every part >= F_min from every part it is not joined to, ports apart (build.py:981-1014; layout.py:1560,1576) -> LayoutCollision refusal before any photomask (build.py:1173-1175); remedies named, no automatic mirror search |
| MV3 _solve_alpha rail too narrow | 6356-6358 | rail cannot be made wide enough for T | REFUSE (GateError) | **EQUIVALENT** | build | layout._root names the failure (layout.py:59-65,211) |
| MV4 _solve_alpha clipping at the floor | 6350-6355,6401-6406 | clipped turns equalised to the common attainable R; shortfall recorded | records (gated by V_turn_closure) | **EQUIVALENT** | build | no clipping: loop widths solved exactly (layout.py:176-216); a loop below the floor shows in side_continuity's loop_narrowest (build.py:140-145), the aspect ratio and the DRC |
| J1 junction edge continuity | 6920-6946 | leg <-> turn cross-section mismatch <= 1e-6 um | ACC-PASS | **EQUIVALENT** | build | one unioned polygon: topology single, side_continuity one strip per leg and side incl. loops, features_found (pipeline.py:625-658) |
| J2 turn rho-closure from the drawn widths | 6953-6983 | abs(R_rail/T - 1) <= 1e-9 and abs(R_in/R_out - 1) <= 1e-9 (unclipped turns) | ACC-PASS | **EQUIVALENT** | build | loop R measured per side vs its target: `validation.loops_measured[].A/B.R_rel`, `fidelity.loop_R_rel_max` (build.py:137-139,1242) reported; their effect gated by Q1 and auditor P3 (<= 0.15) |
| J3 no fluid in a post | 6990-7017 | 0.8 um interior margin | ACC-PASS | **EQUIVALENT** | build | posts are holes of one polygon by construction; holes_ok hard; DRC solid_layer_overlap_with_fluid recorded (drc.py:116-119) |
| J4 minimum corner radius (morphology) | 7024-7029,6857-6882 | open/close at r changes area <= 8e-4 (run only when `rounding` is given; run_acceptance does not) | DEF in acceptance | **EQUIVALENT** | build | rounding_check: min local curvature radius >= process.min_radius_um (build.py:103-123; pipeline.py:652, hard); auditor P10; DRC opening at F_min |
| J5 posts region | 7032-7043 | none in a meander; none above a seam | DEF (never called) | **NOT-APPLICABLE** | - | meander/seam |
| J6 ports and leads | 7046-7065 | exactly INLET, OUT_C, OUT_S+, OUT_S-; 4 leads; each port inside a channel | ACC-PASS | **EQUIVALENT** | build | components_drawn, requested_vs_drawn, port_separation, port_spacing (pipeline.py:633-650, hard) |
| J7 width fidelity | 7068-7074 | abs(w_s drawn - design) <= 1e-9 um on every row (row_geom metadata); w_c non-increasing | ACC-PASS | **EQUIVALENT** | build | measured on the polygon, not metadata: `fidelity.side_eqwidth_dev_max_um` reported, Q1 hard; w_c monotonicity not checked (the tapers are monotone by construction) |
| J8 F from polygons | 7077-7101 | abs(F_em(measured) - F_em(design widths)) <= 2e-3 | ACC-PASS | **PORTED-GATED** | build | Q1 F_abs_dev <= 1e-4 (pipeline.py:661) |
| J9 composition clearance | 7104-7171 | non-adjacent groups >= F_min (distance transform) | ACC-PASS | **PORTED-GATED** | build | part_clearance (MV2), gated before writing |
| J10 cutoff from polygons | 7174-7202 | max abs(delta(f_design, w_c measured)/delta* - 1) <= 0.02 | DEF (never called) | **EQUIVALENT** | build | auditor P4 (kernel at drawn widths) and P5; Q2 on the emergent cutoff |
| J-self self_test n1-n6 | 7225-7311 | shifted rail -> J1, widened rail -> J2, fluid through a post -> J3, missing lead -> J6, floored widths -> J7, butt end -> J4 | DEF (never called) | **EQUIVALENT** | test | test_auditor.py:233-303; test_negative.py; test_units.py:101-107 |
| DXF-D1 ezdxf structural audit | 7414-7420 | 0 audit errors; > 0 polylines | ACC-REP | **EQUIVALENT** | build | P0 (ezdxf read, closed, no bulges, units 13, valid polygon) and the DRC's invalid_rings = 0 (hard); ezdxf's doc.audit() itself is not run |
| DXF-D2 round trip | 7423-7437 | max vertex deviation < 1e-6 um | ACC-REP | **PORTED-GATED** | build | P0 exact equality (auditor.py:126-128,337) |
| DXF-D3 widths from the file | 7457-7467 | all N rows measured; abs(w - design) < 1e-3 um | ACC-REP | **EQUIVALENT** | build | in the build the drawn polygon is measured and P0 proves it is the file; `opencif audit RUN_DIR` measures the file itself (auditor.py:565-617) |
| DXF-D4 F from the file | 7470-7475 | abs(F_em - F_desired) <= 3e-3 | ACC-REP | **PORTED-GATED** | build | Q2 F within F_tol, default 1e-3 (pipeline.py:668-670; spec.py:674-676) |
| DXF-D5 drawn vs design network | 7478-7487 | bulk abs(f ratio - 1) <= 2e-3; abs(F ratio - 1) <= 2e-2 | ACC-REP | **PORTED-GATED** | build | Q1 (every gap 2.5e-3; F 1e-4) |
| DXF-D6 schedule-basis cutoff | 7490-7503 | <= 5e-3 | ACC-REP | **EQUIVALENT** | build | auditor P4 5e-3 vs the kernel at the drawn widths (hard) |
| DXF-D7 adherence and purity from the file | 7506-7529 | reported (always ok) | ACC-REP | **EQUIVALENT** | build | realised_cutoff per gap (build.py:1256); Q2 gates every undeclared gap |
| EMIT emit_device_dxf stamp | 7925-7945 | META text (version, hashes, UTC), .sha256 sidecar | acceptance output | **EQUIVALENT** | build | the record carries the code stamp (pipeline.py:51-68) and the photomask sha256 (auditor.py:549); P0 forbids extra DXF entities |
| OPT1 optimizer criterion: F | 7583-7587,7699 | abs(F_em - F_desired) <= 3e-3 | sweep ok | **PORTED-GATED** | build | Q2 F (default 1e-3) |
| OPT2 optimizer criterion: adherence | 7585-7586,7624-7651,7700 | >= 95 % of filtrate mass from rows with delta_em within +-5 % of delta* | sweep ok | **EQUIVALENT** | build | stricter: Q2 requires every undeclared gap within +-1 % of D on both sides |
| OPT3 sigma_0 compensation | 7605-7621 | fixed point to abs(F_em - F_desired) <= 1e-3; GateError if F_target' <= 0.02 | REFUSE (GateError) | **EQUIVALENT** | build | floor start stops on the delivered F = 1 - P/(1 + 2 S0) (schedule.py:697-700) |
| OPT4 lamina purity | 7654-7672 | 1 - C reported | ACC-REP | **EQUIVALENT** | build | Q2 per gap; realised_cutoff worst_over |
| OPT5 sweep / pick smallest N | 7688-7712 |  | optimizer | **NOT-APPLICABLE** | - | N is free and not minimised (App. B C8 CQ5) |
| FR1 check_fold_routability / G14_fold_routability | 7721-7776 | single track iff Lambda < w_min + 2 F_min; spend of the folds beyond the 8 largest vs band W - 110 um; strict raises | ACC-PASS (FOLD-45 must be refused) | **NOT-APPLICABLE** | - | no fold field is drawn; the designs that would need one are refused in strict mode (straight w_s < w_min, pipeline.py:1053-1063; neg test_negative.py:249-260) |
| V_op Re | 7779-7791,7799 | inlet Re <= 1.0 at 5 uL/min (rho 1000, mu 1e-3) | ACC-REP | **PORTED-REPORTED** | build | gates.P8 at 5 uL/min and the design flow (build.py:911-918,937); auditor P8 on the drawn inlet (auditor.py:505-515) |
| V_op entrance | 7792,7799 | 0.05 Re D_h <= Lambda/2 | ACC-REP | **PORTED-REPORTED** | build | as V_op Re |
| V_op shear | 7793-7794,7800 | gamma_dot = C(alpha) Q/(a b^2) at the inlet <= 1500 1/s at 5 uL/min | ACC-REP | **MISSING** | - | gamma_dot is reported at the design flow only (gates.shear_peak, build.py:837-843); no 5 uL/min value, no 1500 comparison (auditor P8-M14 computes gamma but does not record it, auditor.py:487) |
| V_op pressure | 7795-7800 | trunk-only dP = sum r(w_c) Lambda P Q <= 200 kPa at 5 uL/min | ACC-REP | **MISSING** | - | B3 |
| V_steric D/G | 7805-7811 | D/G < 0.5 | ACC-REP | **PORTED-REPORTED** | build | passage_report: the measured gap opening against 2 D, a warning note when G < 2 D (pipeline.py:521-566,729-733; `validation.passage.classes`); the equality case differs (opencif flags G < 2D, CIF_v2 passes only D/G < 0.5) |
| V_steric D/h | 7805-7811 | D/h < 0.3 | ACC-REP | **MISSING** | - | not computed |
| V_steric delta*/h | 7805-7811 | delta*/h < 0.2 | ACC-REP | **EQUIVALENT** | build | stricter and gated: E11 delta*/h <= 0.10 (A9) |
| V_fab_tolerance | 7814-7867 | 24 draws of N(0, 0.25 um) on every w_c and w_s; yield = share with abs(F - F_desired) <= 3e-3; ok iff yield >= 0.95; adherence mean/min | ACC-REP | **MISSING** | - | no Monte-Carlo and no F yield; the nearest record is robustness_design_network: a uniform bias +-0.5 um and ONE draw of per-segment side noise sigma 0.01 um, reporting cutoff shifts (scenarios.py:32-77; build.py:1268-1280) |
| V_turn_closure | 7870-7881 | zero clipped turns and abs(R/T - 1) <= 1e-9 on both rails | ACC-PASS | **EQUIVALENT** | build | no clipping exists (MV4); measured loop R reported (J2), Q1 gated |
| V_outlet_divider | 7884-7905 | F_div = 2S/(P + 2S) from the march vs the network F_em within 5e-3 | DEF (never called) | **EQUIVALENT** | build | auditor P6b: 2 R_c/(R_s + 2 R_c) at the drawn outlet vs F_drawn within max(8e-3, 4 f_mean) (auditor.py:450-462, hard); P6 abs(F_drawn - F_march) < 0.02; trace eq. 44/75 tests (test_appendix_trace.py:139-150) |
| ACC goldens | 7964-7978 | N = 406/562/731 at w_c0 80/90/100 and abs(F_em - 0.35) <= 3e-3 | ACC-PASS | **NOT-APPLICABLE** | - | CIF_v2's own demo devices (opencif's regressions: test_presets_match_replicas.py, provenance compares) |
| ACC G14 refuses FOLD-45 | 7980-7989 | RoutabilityError expected | ACC-PASS | **NOT-APPLICABLE** | - | FR1 above |
| ACC J-suite on the serpentine-8 fixture | 7991-8007 | J1, J2, J3, J6, J7, J8, J9 and V_turn_closure all ok | ACC-PASS | **EQUIVALENT** | build | the per-build counterparts above, on every design rather than one fixture |

## F. The embedded physics core (_PHYSICS_SRC = ladder_physics.py, CIF_v2 line 3718; line numbers below are in the extracted copy aud/_PHYSICS_SRC.py)

Role "AUDIT" = a check of `audit()`; validators Z1, Z2, T1, W3, W4, Q1, Q5 and R2a require every P check to pass on the demo devices, and compile_unified refuses when P3 or P4 fails (U3).  opencif's port is opencif/checks/auditor.py, run on every photomask (build.py:1332-1335) and made a pass condition by pipeline.assess (pipeline.py:679-697).

| ID / name | CIF_v2 line(s) | what it checks, threshold | role in CIF_v2 | opencif 1.0.0 status | level | evidence (file:line, record field) |
|---|---|---|---|---|---|---|
| F1 physics_march lamina overlap | 333-335 | delta* < w_c/2 | REFUSE (ValueError) | **PORTED-GATED** | build | A5 |
| F2 taper exhausted | 345-347 | w_c > 4 delta* while marching | REFUSE (ValueError) | **NOT-APPLICABLE** | - | canonical (Psi) taper only; opencif's constant-shear taper is bounded below by max(wc_min, 0.3 w_c0) (schedule.py:664-665) |
| F3 march terminates | 328,348-349 | <= 200000 rows | REFUSE (ValueError) | **PORTED-GATED** | build | A11 |
| F4 mass identity in the march | 356-364 | abs(P + 2S - 1) <= 1e-12 | TRIP (AssertionError) | **PORTED-GATED** | test+build | test_appendix_trace.py:93-96 (1e-12); auditor P2 (1e-10) |
| F5a unit width >= w_min | 733-737 | every solved unit width >= w_min | REFUSE (ValueError) | **EQUIVALENT** | build | GS >= w_min (wmin_check) and the drawn head units in the aspect-ratio inventory (auditor.drawn_widths, auditor.py:295-296) |
| F5b wall, slot and pocket rules | 432-434,738-743 | b_wall >= wall_min; Lambda - w >= wall_min; m_off - w/2 >= wall_min | REFUSE (ValueError) | **EQUIVALENT** | build | DRC solid band at F_min (hard); head wall_M defaults to F_min (schedule.py:541-542) |
| F5c corner_r cap | 795-812 | corner_r in (0, 0.45 x min white feature] | REFUSE (ValueError) | **EQUIVALENT** | build | D-section U3 |
| F5d other compiler refusals | 480-482,587-593,662-665,683-686,719-728 | >= 1 grid unit; arm headroom; loop convergence; rail drop > 0; decay law; port width drift <= 5 % | REFUSE (ValueError) | **NOT-APPLICABLE** | - | meander_core compiler internals |
| F6 compiler mirror assertions | 772-793 | derive() reproduces pitch, wall, reaches, widths (1e-6 to 1e-9) | TRIP | **NOT-APPLICABLE** | - | meander_core |
| AUD-P1 kernel anchors | 1036-1041 | B(1) = 0.4217310 +-5e-7; c2(1) = 3.5568 +-5e-4; F_lam(0.3, 1e4) = 0.216 +-5e-6; F_lam(0.5, 1.87) = 0.5 +-1e-12 | AUDIT | **EQUIVALENT** | build | auditor P1 (auditor.py:343-352, hard): B(1) +-5e-7, square duct 28.454154 +-5e-6, g(0.3) = 0.216, g(1/2) = 1/2, x(f) round trip < 1e-12, kappa_g = 4.9138: cubic-pair anchors replace the c2 and F_lam ones |
| AUD-P2 mass identity | 1042-1048 | abs(P + 2S - 1) < 1e-10 on the kernel march at the measured widths | AUDIT | **PORTED-GATED** | build | auditor P2 (auditor.py:373-393): plan identity with S0 < 1e-10, plus KCL, product form, outlet sum, solver agreement |
| AUD-P3 drawn closure (sanity bound) | 1049-1065 | max abs(R_exc/(rho R_c) - 1) <= 0.15 | AUDIT | **PORTED-GATED** | build | auditor P3 (auditor.py:394-418): exact and PE closure, both <= 0.15; parity test_auditor.py:159-163 |
| AUD-P4 network bulk | 1025-1030,1066-1071 | max abs(f_em/f - 1) <= tol_f_bulk 5e-3 over range(1, U - n_term), n_term = max(2, ceil(lambda ln(1/1e-3)) + 1), lambda at q_min | AUDIT | **PORTED-GATED** | build | auditor P4 (auditor.py:419-432); parity test_auditor.py:144-152 |
| AUD-P5 lamina uniformity | 1072-1076 | max delta_em/delta* <= 1.02 | AUDIT | **PORTED-GATED** | build | auditor P5 (auditor.py:433-444), both sides, undeclared gaps; parity test_auditor.py:153-158 |
| AUD-P6 filtrate fraction | 1077-1083 | F_em >= F_target - 1e-3 and abs(F_em - F_march) < 0.02 | AUDIT | **PORTED-GATED** | build | auditor P6 (auditor.py:445-449); parity test_auditor.py:164-166 |
| AUD-P6b corrected closed form | 1084-1092 | abs(F_Psi - F_march) < max(8e-3, 4 f_mean), F_Psi = 1 - Psi(w_N)/Psi(w_0) (witness of the doc-form error) | AUDIT | **EQUIVALENT** | build | auditor P6b (auditor.py:450-462): same tolerance, a different closed form (App. A's outlet identity at the drawn outlet); the Psi form needs the canonical taper, which opencif does not offer |
| AUD-P7 port matching | 1093-1097 | 3x rail run-out mismatch; always ok | AUDIT (report) | **PORTED-REPORTED** | build | auditor P7 (auditor.py:463-471), both side outlets x3; parity test_auditor.py:141-143 |
| AUD-P8 Stokes gates at the operating point | 1098-1111 | inlet Re <= Re_max = 1.0 at Q_in_uLmin (default 5) and 0.05 Re D_h <= Lambda/2; gamma_max printed | AUDIT | **PORTED-REPORTED** | build | auditor P8 (auditor.py:505-515) role report, at the design flow and 5 uL/min; parity test_auditor.py:167-170; gamma_max not recorded; LIMITATIONS 6 explains the departure; CU4 |
| AUD-P9 fabrication floors | 1112-1119 | min unit width >= w_min; boundary wall >= wall_min; min_clear >= wall_min | AUDIT | **EQUIVALENT** | build | auditor P9 (DRC clean at F_min - 0.1, narrowest drawn CIF channel >= F_min, hard) and P9-AR (h/w <= AR_max: hard in strict, reported otherwise) (auditor.py:516-534) |
| AUD-P10 rounded corners | 1120-1148 | corner_r in (0, 0.45 x min(slot, boundary, reserve, min_clear)], >= 0.9 cap in auto; round joins and caps declared; fillet squares = discs = concave corners; 3 port squares | AUDIT | **EQUIVALENT** | build | auditor P10 + rounding_check: the smallest local radius of curvature of the drawn boundary >= process.min_radius_um (auditor.py:535-541; build.py:103-123) |
| F7 audit-driven width correction | 1158-1206 | damped Jacobi on bulk widths, <= 8 iterations, revert if worse | design step | **NOT-APPLICABLE** | - | opencif's closure is exact; its fail-safe margins (scenarios.py:107-155) are a different mechanism |
| F8 run_selftest anchors | 1222-1232 | B(1), 12/B(1), r symmetry, c2(1), c2(inf), g(0.3), F(1/2), inversion at beta 0.53, C(1) | SELF (Z3) | **EQUIVALENT** | test+build | test_kernels.py bit-identity; auditor P1; the c2/F_lam ones N/A |
| F9 run_selftest mutation probes | 1240-1267 | +25 % width -> P3/P4; +60 um reach -> P3/P4; corner_r removed -> P10 | SELF (Z3) | **EQUIVALENT** | test | test_auditor.py:242-294 (gap widened -> P4, P5; side narrowed -> P3; sharp corner -> P10) |

## G. The embedded meander core (_MEANDER_SRC = meander_pattern.py, CIF_v2 line 3716)

Its 19 check_* functions are the Y1-Y19 rows of section D (one validator per check); they are not repeated here.

| ID / name | CIF_v2 line(s) | what it checks, threshold | role in CIF_v2 | opencif 1.0.0 status | level | evidence (file:line, record field) |
|---|---|---|---|---|---|---|
| G-M1 Params/derive refusals | _MEANDER_SRC 272-800 (about 45 ValueErrors) | positive box; >= 1 arm and taper; counts >= 0; grid xor release; fraction ordering; line_w <= min_gap - min_clear; hairpin_w; boundary_w >= min_white; trunk_w >= max(2 min_clear, line_w); profiles end at 1; channel fits its row gap; corner_r <= cap; staircase room | REFUSE (ValueError in build_pattern) | **NOT-APPLICABLE** | - | meander_core parameterisation; opencif's own pre-drawing refusals are spec.py (SCHEMA, positive_values, fillet rule, rectilinear_rows >= 1) and pipeline._prepare's head checks (pipeline.py:1209-1228) |
| G-M2 meander run_selftest and MUTATIONS | _MEANDER_SRC 2177-2524 | 17 mutation probes must each be caught | SELF (U1, U2) | **NOT-APPLICABLE** | - | covered by the U1/U2 rows |

## Every MISSING item

* **B3 total pressure drop** (CIF_v2 1303,1314): dP = sum R_c(i) Q_c(i) along the trunk, reported in kPa.  Now: not computed anywhere (no pressure in any record); the ingredient `design.R_element_SI` (gap 1 -> outlets per unit flow, build.py:1311) is recorded  **Suggestion:** Compute dP = design.R_element_SI x Q_in (gap 1 to the outlets, already recorded; add the inlet stem if wanted) and, on the drawn mask, sum R_c(i) Q_c(i) from the auditor's own ladder solution; report it at the design flow and at 5 uL/min beside 200 kPa, not gated (CIF_v2 never gates it).  Planned in memory for this audit: 7.4-86.9 kPa at the design flows of 7 shipped specs, 0.5-3.4 kPa at 5 uL/min: all inside 200 kPa.
* **B6 F_real = eta F** (CIF_v2 1310): eta = 0.903.  Now: not applied and not reported (documented: LIMITATIONS 5; Pairing.eta kept, kernels.py:157)  **Suggestion:** Either report F_real = 0.903 F_drawn labelled as Dinh's measured outlet ratio (CIF_v2 reports it, never applies it), or keep the documented omission (LIMITATIONS 5); reporting it costs one field.
* **B7 transfer bias (S5, P3)** (CIF_v2 260-269,1306): D_true/D_pred = sqrt(mean c2(beta_cal)/c2(beta(i))) over the rows; "REPORTED, NEVER APPLIED".  Now: not computed; opencif's disposition calls beta_cal "not needed on the cubic pair" (test_appendix_trace.py:344; ALIGNMENT AA12, AG2), which is the glossary's statement about N, while CIF_v2 reports this bias on D for its default cubic pair  **Suggestion:** Needs c2(beta) (the F_lam family opencif did not port): port CIF_v2._c2_exact/c2 (with a bit-identity test as test_kernels does) and report [min, max] sqrt(mean c2(beta_cal)/c2(h/w_c(i))) over the rows, never applied; or write the disposition explicitly (the current text, "not needed on the cubic pair", answers a question about N, not about the bias on D that CIF_v2 reports for its cubic default).
* **B8 sensitivity map (S16-e)** (CIF_v2 1101-1132,1307): dlnF/dln w_s(out), w_c(out), h; dlnD/dln w_c, h, f; dlnN/dln kappa = 2; "MANDATORY REPORT".  Now: not computed (the record has dD_last/dln R_side,out and dln R_ext only: pipeline.py:811,816)  **Suggestion:** Finite differences as CIF_v2 does (d = 1e-4): dlnF/dln w_s(out), w_c(out), h with F = 2/(r(w_s)/r(w_c) + 2) at the last row; dlnD/dln w_c, dln f at row 1 with D = kappa_g x(f) w_c (on the cubic pair dlnD/dln h = 0 exactly, which is itself worth stating); dlnN/dln kappa = 2.  Report-only in validation.gates.
* **B9 wall porosity phi = G/Lambda** (CIF_v2 199-202,1193): reported because the kernels assume a solid wall.  Now: not reported (trivially G/(G + g_p) from the spec)  **Suggestion:** Report G/(G + g_p) in validation.gates (one line).
* **P3 transfer bias reported, never applied** (CIF_v2 1818): band computed and > 0; D uncorrected.  Now: D is uncorrected (true) but the bias is not reported: B7  **Suggestion:** As B7.
* **M17 sensitivity map emitted** (CIF_v2 2048): keys present; abs(dlnF/dln h) < abs(dlnF/dln w_s)/3 with both outlet channels < h.  Now: B8  **Suggestion:** As B8.
* **V_op shear** (CIF_v2 7793-7794,7800): gamma_dot = C(alpha) Q/(a b^2) at the inlet <= 1500 1/s at 5 uL/min.  Now: gamma_dot is reported at the design flow only (gates.shear_peak, build.py:837-843); no 5 uL/min value, no 1500 comparison (auditor P8-M14 computes gamma but does not record it, auditor.py:487)  **Suggestion:** gamma_dot at 5 uL/min = gates.shear_peak.row1_per_s x 5/Q_design (shear is linear in Q), or record the gamma the auditor already computes in P8-M14; report against 1500 1/s at 5 uL/min and at the design flow, not gated: like P8's Re, it fails at every shipped design's own flow (2.3k-18k 1/s on the 7 specs planned here) and passes at 5 uL/min (114-1177 1/s).
* **V_op pressure** (CIF_v2 7795-7800): trunk-only dP = sum r(w_c) Lambda P Q <= 200 kPa at 5 uL/min.  Now: B3  **Suggestion:** As B3 (CIF_v2's estimate is trunk-only at 5 uL/min: sum r(w_c) Lambda P Q).
* **V_steric D/h** (CIF_v2 7805-7811): D/h < 0.3.  Now: not computed  **Suggestion:** D/h = S.KAPPA_G delta*/h, one line in build.gates (or beside passage_report), reported against 0.3; shipped designs <= 0.075 (LIMITATIONS 6), 0.016-0.075 on the 7 specs planned here.
* **V_fab_tolerance** (CIF_v2 7814-7867): 24 draws of N(0, 0.25 um) on every w_c and w_s; yield = share with abs(F - F_desired) <= 3e-3; ok iff yield >= 0.95; adherence mean/min.  Now: no Monte-Carlo and no F yield; the nearest record is robustness_design_network: a uniform bias +-0.5 um and ONE draw of per-segment side noise sigma 0.01 um, reporting cutoff shifts (scenarios.py:32-77; build.py:1268-1280)  **Suggestion:** Re-solve the MEASURED drawn ladder (auditor.solve_ladder on m[Rc, RsA, RsB, RgA, RgB] and the outlet), scaling each resistance by r(w + dw)/r(w) with dw ~ N(0, 0.25 um) per width (as scenarios.py does for its bias), 24 seeded draws; report F mean/std, the yield abs(F - F_target) <= 3e-3 (and <= F_tol), and the share of undeclared gaps within +-1 % (opencif's adherence).  Report-only, as in CIF_v2.  It would be a new model_scenarios entry.

## Items whose status is uncertain

* **B7 transfer bias (S5, P3)** (CIF_v2 260-269,1306), marked MISSING: not computed; opencif's disposition calls beta_cal "not needed on the cubic pair" (test_appendix_trace.py:344; ALIGNMENT AA12, AG2), which is the glossary's statement about N, while CIF_v2 reports this bias on D for its default cubic pair
* **P3 transfer bias reported, never applied** (CIF_v2 1818), marked MISSING: D is uncorrected (true) but the bias is not reported: B7
* **M9 invariance sweeps** (CIF_v2 1914), marked EQUIVALENT: by construction: schedule.build never reads Q_in, mu is the constant MU (schedule.py:36), f is assigned from delta*/w_c; no invariance test
* **M19 Gifford Fig. 2, engine-driven** (CIF_v2 3368), marked EQUIVALENT: partial: the implementers' constant-f march (f 1e-4: P, sigma, F at row 5000; test_appendix_trace.py:76-87) checks the closed form at one row; no i_T or depth-collapse check
* **R1 determinism** (CIF_v2 2462), marked EQUIVALENT: rebuild comparisons in provenance (byte-identical photomasks, CHANGELOG rc8); no determinism test in tests/
* **G24 depth-correct field solve** (CIF_v2 3606), marked EQUIVALENT: no 2-D field solve on the mask; a short across a wall is caught by the hole count (a cut post changes it), the per-gap width measurement and Q1
* **R6b DXF posts = the 2:1 pill** (CIF_v2 5676), marked EQUIVALENT: by construction (posts drawn by stadium_post(g_p, l_g), layout.py:736-749; P0 file = drawing; G measured at mid-wall); l_g itself is not measured (LIMITATIONS 12)

## Where opencif's documentation claims a correspondence the code does not support

| claim | where | what the code does |
|---|---|---|
| The auditor is "not run on the exported mask" | docs/LIMITATIONS.md:50-79 (item 6, "opencif therefore does not run CIF_v2's physics auditor (P1-P10) on its exported mask"); docs/ALIGNMENT.md:91 (AB27); docs/DISSERTATION_NOTES.md:265-278; EVALUATION.md:7,154-158 and docs/EVALUATION.md (same lines); README.md:170-172; foundation/README.md:24-35 (open item D95's rationale) | build.py:1332-1335 runs opencif.checks.auditor.audit on every photomask; pipeline.py:679-697 makes its hard checks a pass condition in every mode and P8-M14 part of the envelope; api.py:82-85 (Result.physics_audit); compose.py:51-53 refuses an element whose photomask failed it.  The record field the docs name, `assessment.status.physics_auditor_P8` "NOT RUN ON THE MASK", no longer exists: the code emits `status.physics_audit` (pipeline.py:574,898,922).  docs/CHANGELOG.md has no 1.0.0 entry (its top is rc8). |
| "a correspondence of purpose, not a port, and the closed-form regression witness (P6b) has none" | docs/LIMITATIONS.md:75-79 | auditor.py ports P2-P8 with CIF_v2's tolerances and test_auditor.py:123-170 checks parity with CIF_v2's own audit() on demos A-C; P6b now exists, with a different closed form (auditor.py:450-462). |
| "every check keeps its purpose, formula and threshold" and the test "requires the same numbers" | opencif/checks/auditor.py:8-11 (module docstring); opencif/__init__.py:19-20 ("checks P1-P10 on it") | Formula changed for P1 (cubic anchors replace c2(1) and F_lam; auditor.py:343-352), P6b (outlet identity instead of the Psi form), P9 (DRC and narrowest drawn channel against F_min, plus P9-AR, instead of unit width >= w_min and wall >= wall_min), P10 (smallest curvature radius of the drawn boundary instead of the declared fillet/cap rendering and the 0.45 x min-white cap); role changed for P8 (an audit pass condition and a compile_unified refusal in CIF_v2, a report in opencif).  test_auditor.py compares numbers only for the network, the terminal window, P3 (its PE part), P4, P5, P6, P7 and P8 (test_auditor.py:1-6,138-170), not P1, P2, P6b, P9, P10. |
| V_steric: opencif computes "neither of the others" (V_steric, V_fab_tolerance) | docs/LIMITATIONS.md:81-87 | Understated for V_steric: D/G < 0.5 is computed as the passage report's gap opening against 2 D (pipeline.py:521-566,729-733), and delta*/h < 0.2 is implied by the E11 gate at 0.10 (build.py:932); only D/h is absent.  V_fab_tolerance is indeed absent.  The paragraph does not mention V_turn_closure, V_outlet_divider, check_fold_routability, the J-suite or the D-suite, whose dispositions appear nowhere in the docs. |
| Every CIF_v2 field "accounted for" reads as full coverage | docs/ALIGNMENT.md:130-157 (sec. 5); test_appendix_trace.py:330-379 | True for the 39 fields of Pairing, Tolerances and Process, but the dispositions do not cover design()'s mandatory reports (dP total, sensitivity map, transfer bias, porosity, F_real: section B, all MISSING without a stated disposition except F_real, LIMITATIONS 5), the 162 validate() validators, or the acceptance validators; this inventory is the first place they are listed. |
| "`cifsched.rr` -> `CIF_v2.r`", "`cifsched.f_of` -> `CIF_v2.g_cubic`", "`measure2.solve3`", "`reference/CIF_v2.py` (byte-frozen kernels)" as the code | docs/ALIGNMENT.md (14 `cifsched` mentions, e.g. AA2, AA3, AA18, AA26); docs/DISSERTATION_NOTES.md:23,29 (9 mentions); docs/LIMITATIONS.md (4); docs/CHANGELOG.md (4) | No module named cifsched or measure2 exists in opencif/; the code is opencif.physics.schedule and opencif.checks.measure, calling opencif.physics.kernels, a port proven bit-identical (test_kernels.py:36-99), and test_kernels.py:99-101 asserts that no package module imports CIF_v2. |
| "opencif's docs do not use spec v5's benign ... O(f^2)" | docs/ALIGNMENT.md:114 (AG9) | opencif/physics/kernels.py:188-190 (entrance_length docstring) uses exactly that argument: "the per-row disturbance is of order f ~ 1e-4, so its effect on a resistance is of order f^2". |
| docs/ARCHITECTURE.md "walks through a build" | opencif/__init__.py:12 | docs/ARCHITECTURE.md does not exist in the tree (not a CIF_v2 correspondence, noted in passing). |

## Method notes

* The numbers quoted for dP, wall shear and D/h come from planning seven specs in memory with opencif's own planner (examples/gifford_physics_linear, wbc_paperclip_teardrop, platelet_serpentine_chip; presets gifford2014_fig4_design, gifford2014_fig6_design, xia2016_rs7_specified, dinh2024_v1_specified; a scratch script of the audit, not kept), not from drawn masks; nothing was written to the repository.
* Each row was checked against the source text of both trees; where opencif's counterpart is a test, the level column says so.
