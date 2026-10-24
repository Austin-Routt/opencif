# Errata for the next unfreeze of the foundation (v114)

App. A, App. B and the glossary are frozen at v114 (`README.md`).  Nothing here is applied to the frozen texts; each
item gives the place, the reason and the proposed text, so that the next unfreeze can apply it in one pass.  The
dissertation's software chapter and opencif's docs use the corrected values meanwhile (`docs/DISSERTATION_NOTES.md`).

## D95 — App. B :4373–4377, the warrant sentence (wording): closed into D106

Recorded in the rc4 round, wording revised in rc5; closed into D106 in opencif 1.0.0, which runs the physics auditor on
every photomask (`README.md`, "Open item D95").  No number moves.

## D106 — name the generator beside the design code (proposed with RC8; revised for opencif 1.0.0)

**Places.**  App. A's title page and table caption (A:50, A:184–187), which name the programs behind the results; App.
B :4373–4377, the software sentence, which speaks of "the targets-to-photomask tool (Prospectus C7)" without naming it.

**Reason.**  The dissertation's software is opencif; the foundation names only the design code `CIF_v2.py`.  opencif
1.0.0 restates the design code's kernels (its `opencif.physics.kernels`, tested to return `CIF_v2.py`'s numbers bit
for bit) and runs the design code's physics auditor, ported, on every photomask; the auditor enforces the design code's
channel and Stokes gates (M14) on the drawn mask and checks and reports, without enforcing, its own stricter inlet limit
(P8: inlet Re ≤ 1, which every lineage device exceeds at its own flow).

**Proposed text.**  A:184–187, after "the design code `CIF_v2.py` (Results 10–14)": "whose kernels the targets-to-
photomask generator opencif 1.0.0 (code stamp `9e4eb55a0affeef3`) reproduces bit for bit".  B:4375–4377: "... warranted by
their acceptance suite and the physics auditor (the design code's, which opencif 1.0.0 runs on every photomask file:
it enforces the design code's channel and Stokes gates on the drawn mask and reports its own stricter inlet limit)
rather than by a defeasible inference."  No claim, label, gate or number moves.

## D96 — the glossary's Gifford row pitch, 46.1 → 45.8 µm (rc6)

**Places.**  Glossary :1865–1866 ("Against the lineage row pitch $\Lambda=G+g_p=46.1$ µm, which the record measures
from Gifford's Fig. 4 ($G=20.5$, $g_p=25.6$)"); the table row at :2536 ("against the Gifford Fig. 4 pitch
$\Lambda=G+g_p=46.1$ µm"); the macro `\LeOverLamGInlet` at glossary :296, `PREAMBLE_FINAL.tex` :296 and
`APPENDIX_B_current.tex` :242 (`\providecommand`).  A source comment at glossary :2469 ("Gifford's MEASURED pitch
46.1") is not printed.

**Reason.**  The 46.1 µm came from Fig. 4's inset read at its 300 µm bar counted at its full-intensity pixels only
(107 px).  The bar is anti-aliased; read with its end pixels counted by coverage it is 107.52 px, and the four post
columns of the inset then give Λ = 45.81–45.84 µm, mean 45.82 µm (`provenance/rc6_scripts/measure_rc6.py`,
`opencif/data/r6/figure_anchors.json`).  With G = 20.5 µm (printed), g_p = 25.3 µm.

**Proposed text.**  :1865–1866: "Against the lineage row pitch $\Lambda=G+g_p=45.8~\si{\micro\meter}$, which the
record measures from Gifford's Fig.~4 ($G=20.5$, $g_p=25.3$)"; :2536: "$\Lambda=G+g_p=45.8~\si{\micro\meter}$";
`\LeOverLamGInlet` 1.31 → **1.32** (\LeInlet / Λ = 60.5 / 45.8 = 1.321) in all three files; `\LeOverLamGOutlet`
stays 1.02 (46.9 / 45.8 = 1.024).  The conclusion (L_e exceeds the row pitch over much of the device; expected
benign; report, do not gate) does not change.

## D97 — the glossary's Dinh row pitch, "near 49 µm, an inference" → 45.1 ± 0.2 µm measured (rc6)

**Places.**  Glossary :1867–1872 ("\emph{Dinh publish no row pitch}: their Fig.~2 gives the gap width, the taper and
the depth but no post width or spacing, so an estimated Dinh pitch near $49$ µm (giving \LeOverLamInlet and
\LeOverLamOutlet) is an inference, not a source figure, and the traceable Gifford value is quoted here instead.") and
:1885–1886, which print the two Dinh-pitch macros (the sentence's text holds); the macros `\LeOverLamInlet` and
`\LeOverLamOutlet` at glossary :295 and :298, `PREAMBLE_FINAL.tex` :295 and :298, and `APPENDIX_B_current.tex` :241
and :244.

**Reason.**  Dinh print no row pitch, but their Fig. 2B shows the posts at a printed 50 µm bar on each panel: six
post-column readings on panels 2–4 give 45.14 / 45.17, 45.09 / 45.09 and 45.11 / 45.05 µm, mean 45.11, SD 0.04 µm,
each panel's bar read with its anti-aliased ends (46.906 px); the tolerance ± 0.2 µm is the bar's reading (± 0.3 %)
plus the readings' spread (`provenance/rc6_scripts/measure_rc6.py`; `opencif/data/r6/figure_anchors.json`,
`dinh2024_fig2b_panel4.post_pitch_um`).  It is a measurement on the source figure (M), not a printed value; the
glossary's 49 µm estimate is superseded.

**Proposed text.**  :1867–1872: "\emph{Dinh publish no row pitch}; measured on their Fig.~2B micrographs (each
panel's own $50~\si{\micro\meter}$ bar) it is $45.1\pm0.2~\si{\micro\meter}$ (giving $\LeOverLamInlet$ and
$\LeOverLamOutlet$), a measurement on the source figure rather than a printed value, so the traceable Gifford value
is quoted here as well."  Macros: `\LeOverLamInlet` 1.23 → **1.34** (60.5 / 45.1 = 1.341), `\LeOverLamOutlet`
0.96 → **1.04** (46.9 / 45.1 = 1.040), in all three files.  The next sentence ("On either pitch the ratio exceeds
unity at the inlet and sits at unity at the outlet") and :1885–1886 hold as written.

## rc7 proposals (from the rc7 foundation check, `review7/foundation_check/FOUNDATION_CHECK.md`, 2026-10-19)

Not applied; recorded so that the next unfreeze can weigh them.

- **D98 — glossary cheat sheet, N (:2507).**  "N: design input" contradicts the spec (S:38 "N is an output", S:325) and
  App. A (A:1310–1311, the march stops at the first row with P(i) ≤ 1 − F; A:1354–1355, "the cutoff fixes δ*, and the
  filtrate target then fixes N"), and the glossary's own subsections (G:1454, G:1741; G:2485 lets the subsection
  govern).  Proposed: "N: design output (the march stops at the target F)".
- **D99 — App. B C7's wording for Dinh (B:993–996, :2974–2975; glossary :1321–1322, :1340–1341).**  C7's band is a
  correct computation at the printed f* (held constant, B:2965–2967, or Xia-scaled, B:2970–2972); "for Dinh's device"
  and "conservative envelope" assume the linear section runs at that f*, which the figures' row count contradicts
  through the closure at the printed widths.  Holding D 2.17 µm with any taper needs ≥ 1,762 linear rows (79.5 mm at
  45.1 µm) at the 1-D outlet identity's F 0.5592; at the measured F 0.505 it needs ≥ 1,494 (67.4 mm), and, with the head
  at f* up to the printed 21 µm seam (App. B's own reading, B:3208–3210), ≥ 80.0 mm in all, against the printed 75.8 mm
  footprint.  Figs. 2A and 3A put the linear section at 514–642 rows: the closure at the printed seam and outlet widths
  with the mass identity (A:1126–1132, A:1141, A:1148) and R13(a)'s invariant lamina (A:1294–1303) turn that into a cutoff
  of 4.28–4.79 µm at F 0.5592, 4.04–4.52 µm at the measured F (opencif register INV-1).  Proposed: qualify the sentence
  as conditional on the linear section running at the printed f*, and name the tension.
- **D100 — App. B :3208–3211, Dinh's head share.**  "which the paper does not print": it follows from the printed 21 µm
  seam width: σ_seam = r(75)/r(21) = 0.0299 of the centre flow per side (0.0282 of the inlet flow per side; the head
  extracts 5.65 % of the inlet flow, 10.1 % of the 1-D F), under the PE closure at the seam's straight segment (A:1148)
  and the 1-D duct law at h 140 µm.  Proposed: give the share, its premises and its source.
- **D101 — spec v5 S:456, the constant-shear ODE.**  Holds only for w_c ≤ h (b = w_c).  For w_c > h the shear
  C(α) Q/(w_c h²) held constant gives d(w_c²)/di = −12 δ*² / (1 − w_c C′/C), against the printed −6 δ*² /
  (1 − (w_c/2C) C′); their ratio is (1 + ε)/(2 + ε), ε = −w_c C′/C ≈ 0.22–0.24 on 300 → 150 µm at 140 µm, so the
  printed form gives about 0.55 × the true rate (confirmed analytically by foundation check 2).  opencif solves the shear
  numerically (cifsched.ShearTable) and is not affected.  Proposed: state the domain, or give the w_c > h form.
- **D102 — App. A :1070–1072, the start.**  A:1071–1072 calls w_s(1) against w_s(0) a row-label difference,
  "immaterial for his constant f_gap".  Gifford's printed start rule chooses a nonzero w_s(1) ("close to or equal to
  zero, depending on the device production capabilities available", p. 4498): that is not a label shift but App. A's
  "other start mode" (A:1118–1119), with σ0 > 0, on which both of opencif's Gifford aligned presets rest.  Proposed: one
  sentence saying so beside the sealed start.

## rc7 proposals from the second foundation check (`review7/foundation_check2/FOUNDATION_CHECK2.md`, 2026-10-19)

- **D103 — glossary cheat sheet, Λ (:2505).**  "G, g_p, Λ & design inputs; Λ = G + g_p" contradicts G:451,
  "Λ(i) ≥ max(G+g_p, 2p [if folded]); equality only for an unfolded tight-packed row" (G:2485 lets the subsection
  govern).  It bears on the question whether a meander lengthens a row: a reader of the cheat sheet would conclude it
  cannot.  Proposed: "Λ ≥ max(G+g_p, 2p [if folded]) (Table A.1); = G+g_p for an unfolded tight-packed row".
- **D104 — glossary :2160–2162 against App. B :3497–3501.**  The glossary calls the single centre-channel inlet (P1.5)
  "a numbered assumption of this work and carries no external anchor"; App. B calls it "a device-level architectural
  boundary condition … not a modelling choice".  (Proposed by the first check, not carried here then.)  Proposed: align
  the glossary with App. B, or qualify App. B as "a premise this work adopts for the published lineage, whose devices
  all show a single centre-channel inlet".
- **D105 — spec v5 § 14.6 (:403–409), L_eq.**  The spec names no reference path for L_eq; against the corner-point
  centreline a square corner's L_eq is negative (opencif rc7: −0.24 to −0.56 w), so "typically 2–6 %" (S:407) and
  "under-resistant" (S:757) are statements about U-turn folds.  "Calibrated once at low Re" (S:409) should allow a
  computed L_eq, marked as a model.  Proposed: one sentence each; optionally define G7's "inner corner" (S:702) as the
  concave fluid corner, a wall's end, as S:371 uses it.
- **Optional, spec v5 G11 (:706).**  "Folded iff w_s(i) < w_min.  No other threshold" makes a printed seam width above
  the floor a departure (Dinh's 21 µm against a 19 µm floor, **A**).  A reproduction note in opencif's docs is enough
  (register INV-1); an erratum only if the spec should allow a declared seam width.

## D107 — App. B C4, F's magnitude: "validated" and "closed by 3-D outlet-bifurcation physics" (opencif 1.0.0)

**Places.**  App. B :1015 ("\textbf{validated} ($F$'s magnitude)"), :1018–1023 ("confirmed to within $1\sigma$ by
COMSOL ... closed by 3-D outlet-bifurcation physics"), and the summaries that repeat the status (:393, :562, :834–840,
:3262, :3324).

**Reason.**  The preprint's 3-D slab lands within 1σ of both of Dinh's means with no parameter fitted, so 3-D outlet
flow is sufficient to account for the 1-D identity's excess; but unstated outlet tubing of 6.5–9.7 cm fits both variants
within 1σ as well (`docs/REVIEW_2026-10-08.md` §3), and so does a side-width fabrication offset s ∈ [5.89 %, 11.86 %]
(the preprint, §IV-C(b)), so the attribution is not unique (`docs/REPRODUCTIONS.md`, Dinh row F).  The slab is also not
Dinh's printed or reconstructed geometry (gap opening 22 µm against the printed 19 µm; wall 19 µm against the 54 µm
measured on Fig. 2B; pitch 49 µm against 45.1 µm; outlets 38 µm after the fifth gap, with no outlet path).  Found by the
independent outline review of 2026-10-23 (B01).

**Proposed text.**  :1015 "\textbf{consistent} ($F$'s magnitude, one device family; the cause of the excess open)";
:1023 "(sufficient to be explained by 3-D outlet-bifurcation flow; unstated tubing or a side-width offset fits as well)".
The status vocabulary (:433) gains "consistent: an empirical magnitude matched by a model within the measurement's
dispersion, without a unique attribution".  No number moves; the separation stays proven.

## D108 — spec v5 §10 and M6: the physical F of a DIRECT (floor) start (opencif 1.0.0; applied as revision 5.1)

**Applied.**  Spec v5 is the project's design document, not a frozen text, so the correction proposed below was made in
place at the author's direction on 2026-10-23 (revision 5.1, the specification's §21, "Fixed in revision 5.1"; md5 in
`README.md`).  Every changed line keeps its place.  The frozen v114 texts need no change: App. A's eq. (67) is the
schedule's bookkeeping.  The design code's docstring (`reference/CIF_v2.py` :18) is frozen and still reads "for ANY
sigma_0".  The record below is kept as it was found.

**Places.**  Spec v5 :292–296 (the corollary "the physical terminal fraction, for any σ₀: F = 1 − P(N)", and "the
other two forms carry the seed"), :672 (M6, "F = 1 − P(N)"), the algorithm's line "F := 1 - P  # PHYSICAL, any
sigma_0" (:574), and the lines that set the stop and the outlet ratio from the same F: :322 (ρ^out = (1 − F_target)/(σ₀
+ F_target/2)), :325 and :569 (stop at the first i with P(i) ≤ 1 − F_target).  The design code's docstring says the same (`reference/CIF_v2.py` :18).  App. A's remark on σ₀
(:1118–1119) and its eq. (67) define the fraction that the schedule removes from the central channel and need no change.

**Reason.**  A DIRECT start draws its side channel at a printable start width with a sealed end, so no side flow
enters it (the spec's own network model: "Q_s(0) = 0 is structural", :424), whereas the march sets every side width
from σ(i) = S(i)/P(i) with S(0) = σ₀ > 0.  The sealed device then carries the schedule's flows scaled by 1/(1 + 2σ₀),
and the entrance rows take up the difference: row 1 alone, off-design by construction (:271), where the pressures
equalize at every row as the closure (§9) assumes, and on the drawn network the first few rows, which opencif's
record declares.  The fraction that the drawn device delivers is
therefore 1 − P(N)/(1 + 2σ₀) = 2/(ρ(N) + 2), the form that the corollary says carries the seed, and not 1 − P(N).
opencif stops a floor start's march on that value, and the drawn network's solve confirms it:
`gifford2014_fig4_as_fabricated_20um` (σ₀ = 2.009 × 10⁻²) delivers 0.90647 against its target 0.90644, where
1 − P(N) would be 0.9027 (`docs/results.csv`; `tests/test_units.py` :65, "floor start stops on the DELIVERED F").
The design code's σ₀ compensation (OPT3) also takes the network's F^em as a DIRECT start's delivered F, by a
fixed-point iteration on its F target.  Found by the independent re-audit of
`docs/ALIGNMENT.md` for opencif 1.0.0 (2026-10-23; rows AA16 and S3).

**Proposed text.**  :292 "**Corollary — the schedule's terminal fraction, for any σ₀:** F_sched = 1 − P(N), the
fraction of Q_c(0) that the schedule removes from the primary"; after :296, "**The delivered fraction.** A FOLD
start (σ₀ = 0) delivers F = F_sched.  A DIRECT start is sealed in the drawn network as well (§15), so it delivers
2/(ρ(N) + 2) = (F_sched + 2σ₀)/(1 + 2σ₀), and a DIRECT design stops its march on that value"; :322 "ρ^out =
2(1 − F_target)/F_target" for the delivered target; :325 and :569 "stop at the first i with P(i) ≤ (1 − F_target)(1 +
2σ₀)"; :574 "F := 2/(rho + 2)  # DELIVERED,
any sigma_0 (= 1 - P when sigma_0 = 0)"; M6 "F_sched = 1 − P(N); if σ₀ = 0, |F − 2S(N)| < ε_id and
|F − 2/(ρ(N) + 2)| < ε_id; if σ₀ > 0, the delivered 2/(ρ(N) + 2) agrees with the network model's F^em".  No shipped
number moves: opencif already stops on the delivered value.


## D109 — the glossary's numerical uncertainty of the 3-D slab: the bound u_num ≤ 5 × 10⁻⁵ is withdrawn (opencif 1.0.0)

**Places.**  Glossary `app:units_validation` (:1098–1109: "Each Dinh variant was solved at two mesh levels a ~5-fold
refinement apart … the reported F is identical at both … That bounds u_num at 5 × 10⁻⁵ … at least \uDoverUnum×
below u_D"); Table A.6, row u_num (:1145, "Bounded, not measured: ≤ 5 × 10⁻⁵"); the macros `\uvalExcessVOne`,
`\uvalExcessVThree`, `\uDoverUnumVOne`, `\uDoverUnumVThree`, `\uDoverUnum` and `\uvalExcessMax` (0.048, 0.011, 660,
720, 660, 0.05) in `PREAMBLE_FINAL.tex` :348–356 (`foundation/` and `reference_glossary/`) and App. B :292–300.

**Reason.**  In `cif-reference-implementation/COMSOL` (tag v3.20), the "extra coarse" models (24,783 and 21,448
elements) store result tables identical to the "coarse" models' to 17 significant figures (Q_c 1.5844077400370972E-10
in both V1 files; 1.3557936033734393E-10 in both V3 files), while their degrees of freedom (138,474 vs 606,563 for V1;
121,788 vs 547,783 for V3), Newton error estimates and stored solutions differ. The tables were not re-evaluated after
the re-solve, so the identity of F says nothing about the mesh; the V3 tables predate both V3 re-meshes (created
2026-04-28). An independent open-source solve made during the reviews (P2-P1, three meshes) gave F 0.9 × 10⁻³ (V1) and
1.6 × 10⁻³ (V3) above COMSOL's 0.5070 and 0.6281; ASME V&V 20 counts such a code-to-code difference as supporting
evidence, not as an estimate of u_num (§2-3.3.4). Found by the reviews of 2026-10-23 and 2026-10-24.

**Proposed text.**  The bound is withdrawn, and no value of u_num is claimed until the coarser models are
re-evaluated. u_num is then V&V 20's two-grid estimate (eq. 2-4-12 with Fs = 3, the factor for two-grid studies and for
unstructured refinement; u_num = GCI/1.15, eq. 2-4-14; the refinement ratio from the element counts, 1.672 for V1 and
1.694 for V3), and the macros are recomputed from it. No conclusion depends on the withdrawn bound: E/u_D is +0.06 (V1)
and −0.58 (V3), and any u_num below 5 × 10⁻³ keeps u_val within 1.2% of u_D, so E/u_val stays +0.06 and −0.58.
