# The texts that opencif cites (frozen at v114, 2026-09-29)

This folder holds the author's texts that the code and the documents of opencif cite. Appendix A gives the step-by-step
derivations of the design equations (Results 1 to 15). Appendix B gives the argument for the design method, with its
thesis, its claims C1 to C9, their premises, and their critical questions. The glossary is the reference for every
variable, unit, and equation, and its preamble `PREAMBLE_FINAL.tex` holds the glossary's numeric values as macros. The
three texts were written as the mathematical foundation of the author's dissertation and were frozen on 2026-09-29
(v114; the freeze record `CIF_ARGUMENTATION_FROZEN_v114.tar.gz` has md5 `c6abcca40fce4dcfdf38e5e6471a5a1a`). Each is
given as its LaTeX source and as the PDF compiled on the day of the freeze. A frozen text is not edited, so the errata
and open items found since are listed in `ERRATA_NEXT_UNFREEZE.md`, with their places, reasons, and proposed texts. The
folder also holds the design specification v5 behind the design code `reference/CIF_v2.py` and the differences from
the previous freeze (v110). `tests/test_appendix_trace.py` asserts the md5s of the four LaTeX files, and it checks that
`reference/glossary_values.json` holds the preamble's own values, using the preamble's copy in `reference_glossary/`.
`docs/ALIGNMENT.md` lists every statement of the three texts that the code implements, relies on, or must not
contradict.

The code and the documents cite App. A by result and equation, numbered as printed in `APPENDIX_A_revised.pdf` (for
example, Result 8, eq. 47). They cite App. B by claim, premise, or critical question (C4, P-C7.1, CQ5), or by line of
`APPENDIX_B_current.tex` (App. B :4373), the glossary by the name of its subsection, and spec v5 by section (§13).

`Ladder_Array_Design_Spec_v5.md` is the project's design specification ("spec v5"), which opencif cites for its exact
closure (§13, equivalence-chain identity), the gap resistance form (§14.4), the segment rule (§1, §14.5), the gates
(§18–19: M13, M14, E11, M16, D4, and the pre-flight γ_GS and w_min = max(F_min, h/AR_max)), f_soft (§19, a warning) and
the supplied constants (§2.3, §20).  Spec v5 is the project's design document, not one of the frozen texts: where it
and the v114 glossary differ, the glossary governs.  Its revision 5.1 (2026-10-23, made at the author's direction)
corrects the terminal fraction of a DIRECT start, which a sealed drawing delivers as 2/(ρ(N) + 2) = 1 − P(N)/(1 +
2σ₀) rather than 1 − P(N): every changed line keeps its place, so the line numbers that opencif's documents cite still
hold, and its §21 ("Fixed in revision 5.1") gives the reason and the changes (`docs/ALIGNMENT.md` S3, erratum D108).  One such place: spec v5 calls L_e > Λ "benign" because the
resistance error is O(f²) (§16-c :463, M15 :681, :795), as does the docstring of CIF_v2's `entrance_length`
("This is BENIGN ... O(f^2) ~ 1e-8"); the glossary says that order "needs the first-order term to cancel, which this
record does not show, so it is not used here as a bound" (the project's open item D52).  opencif
does not use the O(f²) argument: it reports L_e/Λ and gates M13, as both texts say.  The same :795 row rejects the
naive entrance length 0.05 Re D_h as a gate, which is why CIF_v2's physics auditor's L_e criterion is reported, not
gated (docs/LIMITATIONS.md 6).

**Open item D95, closed into D106 (opencif 1.0.0).**  D95 (recorded 2026-10-10 in the rc4 review round, wording
revised in the rc5 round) concerned App. B :4373–4377: the targets-to-photomask tool and the machine-enforced envelope
(Prospectus C7, C8) "are delivered as verified software, warranted by their acceptance suite and the physics auditor
rather than by a defeasible inference".  "The physics auditor" refers back to :4349–4350 ("the physics auditor
separately checks the drawn device's inlet against its own, stricter channel limit"), whose item closes with the tag
"code: CIF\_v2.py reynolds, design (M14); physics auditor (P8); reference implementation, VISCOSITY\_PLASMA"
(:4371–4372).  The auditor is CIF_v2's (its `_PHYSICS_SRC` carries the P8 check).  Up to RC8 opencif did not run it on
its exported mask, so the sentence was true of CIF_v2 and not of opencif, and D95 proposed naming the program: "...
warranted by their acceptance suite (and, in the design code \texttt{CIF\_v2.py}, its physics auditor) ...".
opencif 1.0.0 runs the auditor, ported, on the photomask file of every build (`docs/PHYSICS_AUDIT.md`): its P8 checks
the drawn inlet against the auditor's own stricter limit and records the result, and the design code's channel and
Stokes gates (M14) are enforced on the drawn mask.  The sentence now holds of opencif's auditor as well, with that
difference, which the erratum states rather than leaves to the reader: D95 is closed into **D106, revised for 1.0.0**
(`ERRATA_NEXT_UNFREEZE.md`).  The note on "reference implementation" stands: in the frozen texts that phrase names the
preprint's notebook (App. B :4562; glossary bibliography, "CIF reference implementation v3.20", a computational
notebook; App. B :3492–3493 cites "reference implementation \S 8.7, pe\_consistency\_check"), and the tag's
VISCOSITY\_PLASMA is that program's constant (App. B :4359; `CIF_v2.py` has no VISCOSITY\_PLASMA), so the
frozen texts' "reference implementation" is not `CIF_v2.py`, which App. A calls "the design code". The sentence
carries no argument node (App. B: "The silence is a decision"), so no claim, label, gate or number moves.

**Open items D96 and D97** (rc6, 2026-10-17; for the next unfreeze, not applied): the glossary's two row pitches in
its entrance-length paragraph.  D96: Gifford's Fig. 4 pitch 46.1 → 45.8 µm (:1865–1866, :2536; `\LeOverLamGInlet`
1.31 → 1.32), the inset's 300 µm bar read with its anti-aliased ends (107.52 px, not 107).  D97: Dinh's pitch, "near
49 µm, an inference" → 45.1 ± 0.2 µm measured on Fig. 2B (:1867–1872, :1885–1886; `\LeOverLamInlet` 1.23 → 1.34,
`\LeOverLamOutlet` 0.96 → 1.04).  The paragraph's conclusion holds on either value.  Places, reasons and proposed
texts: `ERRATA_NEXT_UNFREEZE.md` (which also points to D95).

**Open items D98–D105** (rc7, 2026-10-19; for the next unfreeze, not applied), from two independent checks of rc7's
decisions against these texts: the cheat sheet's N (D98) and Λ (D103) against their subsections; App. B C7's wording
for Dinh (D99); Dinh's head share, which the printed seam width fixes (D100); spec v5's constant-shear ODE, valid only
for w_c ≤ h (D101); Gifford's printed start rule beside App. A's two start modes (D102); the glossary's P1.5 against
App. B (D104); spec v5 § 14.6's L_eq, its reference path and sign (D105).  Places and proposed texts:
`ERRATA_NEXT_UNFREEZE.md`.

| file | md5 |
|---|---|
| `APPENDIX_A_revised.pdf` | `a972b6e3b1c595cb4ffffd175f542d13` |
| `APPENDIX_A_revised.tex` | `ab28dc46e4e184bc1b0c6d689a1ee129` |
| `APPENDIX_B_current.pdf` | `af17d9db8481b17ba2fcfe3f70a70e4f` |
| `APPENDIX_B_current.tex` | `69b565d0cb2ef73a12778cb2df677732` |
| `PREAMBLE_FINAL.tex` | `3c075af8f63bb1d735f4b43c64f32d14` (unchanged since v108) |
| `glossary_official.pdf` | `9b656101ac6f00b14d2e57164a18e0b5` |
| `glossary_official.tex` | `e9dc24fdead748e9df9e9eebd608a62e` (v114, which reads "Depth factor, CFD" for K_depth^CFD, v114's D90; an earlier copy reads "measured" there) |
| `Ladder_Array_Design_Spec_v5.md` | `8cc1cba8979e6e3a200bc867f99a2039` (revision 5.1; revision 5 was `89ecfb3730bbaac632dff60b2483b5bb`) |
| `v110_to_v114_APPENDIX_A.diff` | `49b6fa487cc52533e788590122c8365f` |
| `v110_to_v114_APPENDIX_B.diff` | `1c4289ea0d8f17f802fdaacddb90f07c` |
| `v110_to_v114_glossary.diff` | `e610791c85e4b810aba88fa15c972162` |

The design code that App. A names, `reference/CIF_v2.py` (v2.1.0), is the project's `CIF_v2.py` byte for byte (md5
`0e0961a50441afb94671b4470d41bace`, asserted by `tests/test_kernels.py` and `tests/test_appendix_trace.py`). The package
does not import it, and the tests compare the package with it.

## v110 → v114

rc1 (code c15f4de265871445) was audited against v110 (A `327aaa73`, B `4279dd70`, glossary `12be7f61`).  The three
`v110_to_v114_*.diff` files are the full differences (changed lines: A 20, B 175, glossary 192).  They were read in full:

* App. A equation labels and numbering are unchanged (88 display environments, identical `\label` set; equations
  10, 36, 44, 67, 69, 75 and 76 print identically), so every equation number ALIGNMENT cites holds.  The changes are
  wording: the O(f²) remark withdrawn, "Fig. 6 platelet-enrichment device", a symbol note, the residual "is attributed
  to" (was "lives") the three model gaps, and w_ch for the channel width in R14.
* App. B: the thesis's envelope sentence now reads "the design tool checks: it refuses a design that leaves the plateau
  or the Reynolds gates, and it reports, without refusing, a design over its length budget" (v110 said "checks and
  enforces"; its Claim field and CQ5 already said the rest).  C8 now attributes the zero ∂D/∂N to the plateau
  schedule (P-C8.1, S1), with C4 supplying that the two targets answer to separate laws; Leg 4 adds "This appendix
  reports no device fabricated to such a schedule".  Other changes reframe or tighten arguments without touching
  anything the code uses: the β law as thin-lamina versus linear (was parabolic versus Hele-Shaw), the depth factor
  divided out as the CFD value, 0.3516 explained as the product N·f_gap, the crowding argument weakened, "within the
  model's premises", "numerical sweep band", Inglis's percentages.
* Glossary: the design code "carries both pairs: its `design` function defaults to the cubic pair" (AG2); the M14
  sentence on the gap's own Reynolds number, (w_c + h)/(G + h) above the flux product (AG10); "Gate instead on the
  perturbation, which the design specification identifies with the design code's gap-dominance gate (M13)" (AG9); the
  O(f²) bound withdrawn (opencif never used it); η listed as measured (as opencif's docs say); the
  `EXTERNAL` marker defined.

No changed statement contradicts the code; three rows gained direct support (AG2, AG9, AG10).  The thesis sentence on
the length budget was already in v110's Claim field and CQ5; rc1 did not obey it (its chip-area check failed a design),
which rc2 corrects (ALIGNMENT AB18, AB26).
