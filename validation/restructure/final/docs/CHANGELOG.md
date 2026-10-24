# Changelog

Earlier rounds (rc1–rc5) are recorded in `provenance/RC*_VERIFICATION*.md` and `provenance/README.md`; the review
history is `docs/REVIEW_2026-10-08.md`.

## 1.0.0 (2026-10-23), code `{{STAMP}}`, tests `{{TESTS_STAMP}}`

Release 1.0.0 turns RC8 into a package that an undergraduate can read and that the dissertation can cite, and it runs
the design code's physics auditor on every photomask. It changes no design and no photomask, because every shipped
specification plans the same row schedule as in RC8 and every photomask is RC8's, byte for byte ({{MASK_RESULT_SHORT}};
`docs/VALIDATION.md`). RC8's records, traces, figures, and register remain the dissertation's evidence for what they
show, and 1.0.0 reproduces RC8's photomasks and adds the audit.

1. The package is restructured (`docs/ARCHITECTURE.md`) into `opencif/physics/` (`kernels`; `schedule`, which is RC8's
   `engine/cifsched`; `network`; and `scenarios`, which is `robust`), `opencif/geometry/` (`shapes`, which is
   `geometry`; `layout`; `depicted`; `dinh_inlet`; and `dxf`, which is `io_dxf`), `opencif/checks/` (`measure`, which is
   `measure2`; `drc`, which is `drc_audit`; `auditor`; `ruler`; and `design_code`), `opencif/figures/` (`proof`, which
   is `renders`, and `drawing`, which is `render`), and `opencif/published.py` (which is `devspecs`). The package uses
   package imports throughout and no longer manipulates `sys.path`. The restructure is reproducible, because the
   scripts in `validation/restructure/` (`make_tree.sh`, with `patch_final.py` for the changes made after the voice
   pass) turn RC8's frozen tree into this one with every edit an exact match.
2. The kernels are restated (`opencif/physics/kernels.py`). The design code's duct law, skim cubic and its inverse,
   laminar coefficient c₂, shear coefficient, entrance length, stadium post, and supplied constants (Pairing,
   Tolerances, and Process) are written to be read, and `tests/test_kernels.py` proves them equal to those of
   `reference/CIF_v2.py`, bit for bit, over dense grids. The package no longer imports `CIF_v2.py`, which moves to
   `reference/`, byte-frozen (md5 `0e0961a5…`), for the tests.
3. Every photomask receives the physics audit (`opencif/checks/auditor.py` and `docs/PHYSICS_AUDIT.md`). The design
   code's auditor (P1 to P10) is ported and run on the photomask file read back from the DXF, together with P0 (the
   file is the drawing), R0 (the auditor's own ruler, `opencif/checks/ruler.py`, measures every ordinary straight row
   again from the file and holds the generator's measurement to it within 10⁻³), R6b (every post of the straight rows
   is the process's post, which applies the design code's acceptance check to every post), and M14 (the design code's
   Reynolds and Stokes gates on the drawn mask). A failed hard check fails the build (exit 2), M14 belongs to the
   validity envelope (exit 3), and P7 and P8 (CIF_v2's stricter inlet limit) are reported. P4 compares the drawn
   fractions with the intended network, which is the kernel for a design without declared departures, because a
   declared departure such as Xia's loop pull-in shifts the fractions of the rows downstream of it without decaying;
   the deviation from the kernel is reported beside it. The command
   `python -m opencif audit RUN_DIR` repeats the audit from scratch. The record's field
   `assessment.status.physics_auditor_P8` ("NOT RUN ON THE MASK") becomes `assessment.status.physics_audit`. On the
   shipped records the audit changes one status and no exit code: `m14_paperclip_const_nomeander_floorstart19`, which
   already fails Q2 (exit 2, as in RC8), is now also outside the validity envelope, because its undeclared floor-start
   rows exceed the Stokes gate on the drawn mask (Re_gap 0.098 and Re_delta 0.18, against 0.01), and its verdict names
   the audit's failure (P5). Every other record keeps RC8's status (`validation/RELEASE.json`, records).
4. The design code's reports are added (`opencif/checks/design_code.py`). The mandatory reports of its `design()`
   (pressure drop, F_real = ηF, transfer bias, sensitivity map, and porosity) and the operating validators of its
   acceptance run (`V_operating_regime` at 5 µL/min and at the design flow, `V_steric`, and `V_fab_tolerance`) are
   computed by its own formulas and written to `validation.design_code_reports`, and, as in the design code, none of
   them is a pass condition. `docs/CIF_V2_REQUIREMENTS.md` lists every requirement of `CIF_v2.py` (282 rows) with its
   status.
5. A teaching module is added (`opencif/teaching.py`). It follows the whole path for a straight ladder in fewer than
   300 lines, from the design to the audited photomask, and its design equals the generator's plan number for number
   (`tests/test_teaching.py`). The tutorial is in `notebooks/`.
6. The tests gain `test_kernels`, `test_auditor` (the port against CIF_v2's own auditor on its demonstration devices,
   and fourteen seeded faults, which exercise every hard and envelope check except P2 and P6b), `test_teaching`, and `test_design_code` (each report against CIF_v2's own
   function, with the design code's determinism R1 and invariance M9 as tests). In `test_appendix_trace`, η and β_cal
   are now reported fields (each asserted equal to CIF_v2's), the disposition of L_eq names both of its replacements,
   and the glossary's numeric macros come from `reference/glossary_values.json` when the dissertation's texts are not
   present. In `test_composition`, a case that declares a FAIL on Q1 or Q2 may now also fail the physics audit, but
   only on the audit's checks that ask the same question (P4 for Q1; P5, P6 and P6b for Q2) and never on its file,
   code and geometry checks; the suite's floor start (case 00) fails the audit's P5, as m14 does. The suite predated
   the audit and was first run on this release during its validation. The build records carry the tests stamp with
   which they were built (`5d83c0644e5bb4e3`), from before this change.
7. The stamps hash every file's path with '/' and in the same order on every system. The code stamp's value on Linux
   is unchanged by this rule, whereas the tests' hash was system-dependent in RC8.
8. The comments and docstrings of every module and test are rewritten in the author's voice, following the style sheet
   derived from Routt et al. (2023). The abstract syntax tree of every file, with its docstrings removed, is unchanged
   by the rewrite, and every code token keeps its text and its column (`validation/restructure/voice_check.py`). Review and check identifiers in the comments are kept, and
   release tags that only dated a change are dropped.
9. The documents are revised. There is a new `README.md`, together with `docs/ARCHITECTURE.md`,
   `docs/PHYSICS_AUDIT.md`, `docs/VALIDATION.md`, and `docs/CIF_V2_REQUIREMENTS.md`; RC8's README becomes
   `docs/REFERENCE.md`; and `docs/EVALUATION.md` is rewritten for 1.0.0 as a single copy. Every statement that the
   auditor is not run on the mask is corrected (LIMITATIONS 6, ALIGNMENT AB27, DISSERTATION_NOTES section 4, the
   reference README, and foundation/README's D95), as are the statement that the auditor's other checks have only
   counterparts of purpose, the statement that V_steric and V_fab_tolerance are not computed (LIMITATIONS 6), RC8's
   module names in the documents that describe the current code (the history documents keep their releases' names),
   and a kernel docstring's O(f²) argument, which ALIGNMENT AG9 says opencif does not use. D95 is closed into D106,
   revised for 1.0.0 (the dissertation's errata list; DISSERTATION_NOTES section 4), and D107 records that App. B's
   "validated" for F's magnitude overstates a comparison whose cause is not unique (`docs/REPRODUCTIONS.md`).
10. The repository ships the texts that the code cites: Appendices A and B and the glossary (frozen on 2026-09-29, as
    LaTeX sources and as the PDFs compiled at the freeze), the glossary's preamble, the design specification v5, the
    errata list, and the differences from the previous freeze, in `foundation/` (`foundation/README.md`), so that
    `tests/test_appendix_trace.py` checks the texts themselves in every copy of the repository.
11. A plan that cannot be completed is a spec error (exit 6) with the design's numbers, rather than a traceback
    (exit 1). RC8 converted a closure without a solution in this way, whereas a march that needed more than
    400,000 rows (for example, D = 0.25 µm on a 100 µm center channel) stopped with the march's RuntimeError.
    `tests/test_negative.py` now seeds that case. The figures' titles separate the name from the targets with a
    colon. Neither change touches a shipped design or photomask.
12. The auditor's ruler finds the trunk's corner at each gap from the posts' tips (`ruler.wall_corners`) instead of
    taking the mean of the two neighboring tips, which missed a change of slope by a quarter of it. The rebuild of
    the release found the difference: where a shear taper meets its floor, the mean placed the corner about 4 nm off,
    and the auditor's network then put platelet_serpentine_chip 0.569% from the intended network (P4, tolerance
    0.5%), although the generator's own measurement of the same photomask agreed with it to 8 × 10⁻⁵. With the
    corners, that design's largest bulk deviation is 0.367%, at a U-turn, where the generator's network shows the
    same. The change is in the auditor only, so no design and no photomask changes; `tests/test_auditor.py` adds a
    synthetic taper with a change of slope (part 4), and `docs/PHYSICS_AUDIT.md` explains the construction.
13. Documents only. `docs/DISSERTATION_NOTES.md` item 6 and `docs/LIMITATIONS.md` item 24 quoted rc6's trace medians
    as the traces of the shipped masks; they now point to the shipped masks' medians, which
    `provenance/rc7_traces/README.md` tabulates (Gifford's Fig. 4 inset, leg 4, at 0.61 µm on the nested preset rather
    than rc6's 1.21 µm). The envelope summary of `docs/DISSERTATION_NOTES.md` section 4 now gives Fig. 5C and 5D's
    gap dominance as their records state it (M13 0.083 and 0.138, not 0.080 and 0.133) and separates the Fig. 4
    design preset (M13 0.031) from the aligned one (0.048). The JSON configurations of `provenance/rc6_traces/` name
    their inputs by paths relative to the rc6 working folder rather than by the build machine's absolute paths.
14. Documents, tests' workflow and provenance only, after two independent reviews of the 1.0.0 code (2026-10-23).
    The documents say what the code does where they said more: the width check takes the head units at their planned
    widths (`docs/LIMITATIONS.md` item 7); strict mode gives exit 3 when an envelope check fails only on the drawn
    mask (`docs/PHYSICS_AUDIT.md`, `docs/ARCHITECTURE.md`, and the options sheet); a head unit's measurement sees a
    local narrowing weakly (item 12); and fourteen seeded faults exercise every hard and envelope check except P2 and
    P6b. The findings that need a code change are listed in item 31 for the next release. On Windows, 1.0.0 runs in
    Python's UTF-8 mode (README, "Installation"), and the workflow runs the suites so. `provenance/rc4_scripts/` and
    `rc5_scripts/`, which the documents cite, are distributed again, with the release candidates' folders written
    relative to the working folder, and `provenance/README.md` names what is not distributed. `validation/` gains
    RC8's photomask hashes, both releases' plans, and the record comparison (`rc8_comparison/`), and `NOTICE.md` states
    what the MIT License covers. The tutorial's window of section 5 keeps only the polygons of the drawing, and its
    section 7.1 compares the two audits' verdicts and lists the values that differ. The device kits found one more
    restructure error outside the package: the import rewrite had turned five trace scripts' tuples of folders into a
    single name (`for p in (R):`), so `run_xia.py` and `run_xia_more.py` could not import opencif in the Xia kit; the
    five scripts of `provenance/rc7_traces/code/` and their copies in `rc6_traces/code/` read `(R,)` again, and the kits'
    manifests were rewritten. An independent re-audit of every row of `docs/ALIGNMENT.md` against `foundation/` and
    the 1.0.0 code (criterion 5 of `docs/EVALUATION.md` section 8) relabelled AA16, whose floor start departed from spec
    v5 §10 and M6 (new row S3, erratum D108), and made ten smaller corrections in the matrix. At the author's
    direction, the specification itself was then corrected: its revision 5.1 states the delivered fraction of a DIRECT
    start and stops on it, with every changed line in its place (`foundation/README.md`). Running `docs/EVALUATION.md`
    sections 4 and 5 as printed (criterion 6) found section 5's comments on two refusals incomplete: they now also
    name the refusal on w_min at row 1 (since rc4) that comes before the plateau and Reynolds refusals. `CITATION.cff`,
    `pyproject.toml` and README name the repository, https://github.com/Austin-Routt/opencif. The six publishers'
    figures that the trace suites and the device kits read are distributed in `notebooks/figures/`, with the Royal
    Society of Chemistry's permission and under CC BY 4.0 (`NOTICE.md`); the kits' source notes say so, and their
    manifests were rewritten after the kits' runs (`docs/VALIDATION.md`, section 6).

Some things are unchanged by design, so that every photomask stays RC8's. They are the physics, the schedule, the
closure, the layout, and the drawing, together with three known behaviors left for a later release (LIMITATIONS 30).
A head-less floor start outside reproduce mode is noted at the plan stage and fails Q2 at the build, an F overshoot
beyond the tolerance is noted and fails Q2 unless `physics.trim_F` is set, and the leg fitter's collision search can be
slow. The British spellings of twelve specification and record keys (for example `centre_channel` and
`realised_cutoff`) and of the program's printed messages are also kept, because renaming the keys would break every
preset and record and the messages are stored in the records and compared by the tests; accepting American aliases on
input is left for a later release.

Validation: {{VALIDATION_SHORT}}
## 1.0.0rc8 (freeze 2026-10-21), code `593d977fc83ec217`, tests `984445646aa4aa85`

RC8 releases what the adjudication of the two third-party reviews of rc7 (`RC7_REVIEWS_ADJUDICATION.md`, 2026-10-20)
recommended, with the decisions D1–D5 taken as recommended and D6 from review 3 of the RC8 candidate
(`RC8_CANDIDATE_REVIEW.md`, adjudicated in `RC8_FINAL_REPORT.md`, 2026-10-21). The physics, the schedule, the closure,
the layout and the drawing are unchanged; the foundation (App. A, App. B, glossary v114, spec v5) is not edited.

1. **The Gifford Fig. 4 reproduction is the nested reading** (D1, D6; register GIF4-1, INV-2). `gifford2014_fig4_aligned`
   keeps its name, every measured value (G 17.2 / g_p 28.6 µm, the legs 818 µm apart, the registration) and the 8 µm
   start (an inference from Gifford's 10 : 1 limit), and takes the nested layout of the caption ("the second and fourth legs"), the schematic and Fig. 6,
   which the inset's arrows suggest (the weakest line: they may both mark leg 4). Reading B (rc7: legs in path order) is rejected by the inset's own taper direction
   (V13) and kept as the variant `provenance/rc7_traces/variants/gifford2014_fig4_pathorder.yaml`. Not reproduced, and
   stated: leg 2's side channels, 17.5–18.1 µm per side wider on the image than the printed schedule draws them.
2. **Input validation and the two rc7-final defects** (D2: an RC8, not a documents-only rc7.1): the candidate patch
   reviewed three times (the independent patch review, review 3, this release's code review), plus review 3's E9 and
   the code review's B1 and n1, and its re-check's N1 (NumPy's bool).
3. **The traces carry A7** (coverage and censored null ratios in every record and caption; claims judged on them), Xia
   Fig. 3d (both masks) and 3g (the depicted mask) move from claimed to recorded (D5), and V13, the taper-direction
   sign test, is part of the trace code and of `tests/test_trace_anchors.py`.

Every record was rebuilt at the RC8 code and compared with rc7 final (`provenance/rc8_scripts/compare_rc7_rc8.py` →
`logs/compare_rc7_rc8.{md,json}`): of the 53 records (51 builds, one of them refused as in rc7, and the 2 arrays), 52 are rc7's in every number apart from what RC8 adds (`validation.drc.file_final` on the 50 full builds; the K13 note on `m14…_floorstart19`; `checks.inter_copy_port_spacing` on the 2 arrays), with 51 of the 52 photomasks byte-identical and none renamed; the one that changed is the Fig. 4 reproduction, with its input (D1: photomask `edaf92626d574aea…`, N 3543 as before).

All suites pass at code `593d977fc83ec217`, tests `984445646aa4aa85` (`logs/test_*.log`): `test_appendix_trace` 71, `test_options_sheet` 20, `test_bend_fd` 3, `test_units` 32, `test_lattice` 42, `test_presets_match_replicas` (5 presets), `test_negative` 63, `test_composition` (22 cases and the extras, in two parts), `test_rc7` 112, `test_rc8` 67, `test_figure_anchors` 42, `test_trace_anchors` 339 (with the rasters; "339 checks, 0 failed"), every one passing (exit 0); from the evaluation archive `test_rc8` reports 63 checks ok and 2 not run (the A group's two arrays) (INCOMPLETE, exit 4, by design) and `test_figure_anchors` 25 ok and 17 not run (24 and 17 without the rasters).

### Code

| file | change | closes |
|---|---|---|
| `opencif/spec.py` | D must be finite. `positive_values`: 37 lengths, depths, flows, radii, tolerances and counts (nested keys included: the outlet tubing's bores, the ring's radii) must be finite and > 0, and nine may be 0 (the fillet, the gap-entrance coefficient, the chip margin, the three corner and junction radii, the inlet port gap, the flare start, the tubing length); a true/false value on any of them is refused (review 3, E9); each refusal names its field. Newly refused, accepted by rc7: a negative or zero flow, an array count of 0 or below, a negative chip; `chip: (0, 0)` (rc7: no chip area) and `array.bus_width_um: 0` (rc7: the default bus) (patch review n1); a negative or zero depth, which rc7 refused only from inside the planner without naming the field, is now refused by name. Two conflict pairs: `meander.enabled` against a contrary `reproduction.head` (both ways), and an explicit width that sets the head's arms (`narrow_um` for the narrow style, `channel_um` for pack and arms) beside a different `reproduction.w_arm` (D3, a policy: one source for the head's widths), one message per conflicting head | K1–K5, K7 |
| `opencif/api.py` | every numeric argument (D, F, the ratio, the chip, the taper widths, the flow) converted through `SpecError`; an unknown process and `seal` without a head are `SpecError` (rc7: `ValueError`); `array=True` is refused (rc7: one copy; review 3, E9) | K4, K6 |
| `opencif/pipeline.py` | `plateau_preflight`: in strict mode a headed design below the κ plateau is refused before the head's sealing search (Review 1's s04 design: the frozen rc7 code gave no answer within 1200 s; the candidate refused it in 49–52 s), its refusal carrying the plateau's advice and a `gates_note`; `name_mask` names the file just assessed and removes the other verdict's stale file; `validation.drc.file_final` beside `validation.drc.file`; a floor start outside reproduce mode notes before drawing that its entrance rows are not declared and, where they over-extract, FAIL Q2 (K13), and the sealed-start refusal's remedy says so; a plan whose F overshoots the target beyond `checks.F_tolerance` notes before drawing that Q2's F check is expected to FAIL, naming `physics.trim_F` (K14) | K8, K9, K12, K13, K14 |
| `opencif/build.py` | each outlet-stretch extension try restarts from the spec's `cap_frac` and side inlet width; the kept try's values are restored with its stretch; a later try that cannot meet the floor ends the loop as an outlet that no longer closes does | K10 |
| `opencif/compose.py` | clears the array's own six output names before composing; `inter_copy_port_spacing()`: the 2 mm rule from each copy's own inlet port as drawn (or the shared port) to every other copy's outlet ports, gated with the other pre-write checks; the neighbouring copies' outlet ports reported (D4) | K8, K11 |
| `opencif/api.py` (RC8 code review, B1, N1) | `generate()`'s numeric arguments refuse a true/false value, Python's or NumPy's (any value whose dtype kind is `b`; NumPy 2 names its scalar type `bool`, so a test of the type's name missed it: N1) (rc7 and the candidate turned `flow_uL_min=True` into a flow of 1.0) | E9 completed |
| `opencif/__main__.py` (RC8 code review, n1) | the command line's `--array 0` is passed on and refused, as `array=0` is (rc7 dropped it as "no array") | |
| `opencif/__init__.py` | `1.0.0rc8` | |

Not changed, deliberately: the fitter's early stop (Review 1 of rc7; it refuses a valid design rc7 fits), process
isolation with hard timeouts and a decoupled stamp check (Review 2 of rc7), the DRC's audit sidecar, and the physics.

### Tests (tests hash `984445646aa4aa85`, 12 files)

- `tests/test_rc8.py` (new): 67 checks in groups V (input validation, the command line's `--array 0` included), C (conflicts), P (the plateau
  pre-flight), S and Q (the plan notes K13, K14), N (`name_mask`), I (the inter-copy rule), L (a fit), D (`drc.file_final`
  on a full build) and A (both shipped arrays recomposed over planted stale files); each fix with a control that rc7's
  behaviour on valid input is unchanged. Without the run pickles in `runs/` (the evaluation archive) its A group cannot
  run: INCOMPLETE, exit 4.
- `tests/test_rc7.py`: GIF4-1's checks follow D1/D6 (the aligned preset is the nested reading: leg 4 ending 2 mm short,
  the 8 µm start, the inset's gap; reading B is read from its variant file).
- `tests/test_trace_anchors.py`: A7 (the censored null ratios, coverage reproduced to ±0.002 and, for a claim, ≥ 0.5)
  and a sixth layer, V13 re-run from the rasters, the trace records and the DXFs, every verdict, position count and
  slope held to its anchor.

### Presets and variants

| item | change |
|---|---|
| `presets/gifford2014_fig4_aligned.yaml` | the nested reading (D1): no `layout.order` (nested), `reproduction.pc_last_short` 2000 (leg 4 ends 2 mm short of the others, as the inset's leads place it), the 8 µm floor start, G 17.2 / g_p 28.6 µm; built: exit 0, REPRODUCED with declared departures (floor start), N 3543, F 0.90650, 99.5 % of F within 1 % of D, photomask `edaf92626d57…` (as the candidate's and review 3's). Its header gives the four lines of evidence and, after them, what it does not reproduce (review 3, E4) |
| `provenance/rc7_traces/variants/gifford2014_fig4_pathorder.yaml` (new) | rc7's `gifford2014_fig4_aligned` (reading B), renamed: a rejected reading, built, traced and shown, not a reproduction |
| `provenance/rc7_traces/variants/gifford2014_fig4_start51.yaml`, `…_fold20.yaml` | header comments only (fold20 is a path-order figure; fold20's `paper.assumed` also quoted, the code review's m6, which changes its record's `input` and no engine value); both rebuilt at the RC8 code, byte-identical to rc7's (start51 `3c72b2d7…`, fold20 `b3fcda80…`; the path-order variant `e8f68a39…` is rc7's aligned photomask byte for byte) |
| the two presets that quote "10 : 1" in `paper.assumed` | the text quoted, so that YAML reads it as text, not as a mapping (the aligned preset and the path-order variant; no other spec has it) |

### Traces (`provenance/rc7_traces/`; not hashed)

- A7 (Review 1 of rc7's patch, verified): every record carries `metric.coverage` (found / eligible samples),
  `metric.n_eligible` and the censored null ratios; every caption states them with the chamfer null; claims are judged
  on the censored ratios and need coverage ≥ 0.5. D5: Xia Fig. 3d (both masks) and 3g (the depicted mask) are recorded,
  not claimed (coverage 0.42, 0.43 and 0.48; censored null ratios 1.00, 1.00 and 0.94).
- V13 (`code/taper_sign.py`, written for opencif): the sign of the image's outer wall-to-wall width slope against the
  mask's along a leg, with an eligibility/coverage rule (the wall-search window must exceed the largest wall offset,
  review 3's E6: ±30 µm on the neighbour, the code review's m3) in place of a fixed position count, which would make
  leg 4's valid control (about 80 positions) inconclusive; a standard error that is the larger of a Newey-West one (Bartlett, lag 100 µm) and the
  independent-sample one (the residuals of neighbouring positions are correlated: the code review, M1, N2); the DXF checked against the traced mask's sha256 (m4); no verdict
  on a mask that does not taper in the window (n2). Required verdicts (`make_anchors.py`): the aligned (nested) reading's neighbour PASS, reading B's
  neighbour FAIL, every leg-4 control PASS. On the inset: the nested reading's neighbour (leg 2) PASS, −0.00240 ± 0.00059 µm/µm on the image against the mask's −0.00257 (4.1 SE); reading B's neighbour (leg 3) FAIL, −0.00254 ± 0.00062 against +0.00301 (4.1 SE); leg 4's control PASS on all five masks (11.6–11.8 SE). Each SE is the larger of the Newey-West one (lag 100 µm) and the independent-sample one.
- The Fig. 4 trace's caption reads the neighbouring leg's walls from the record (review 3, E13).
- `run_all_rc7.sh` builds and traces the path-order variant and runs V13; `make_anchors.py` writes the V13 block.

### Documents

- Register (`RC8_DISCREPANCIES.md`, rc7's with the RC8 entries): GIF4-1 and INV-2 (the decision, its evidence, V13,
  leg 2's open width with its size: review 3, E3), XIA-1 (Xia's own words on the taper's method), XIA-9 (new: A7, D5),
  BUILD-2 (K7, K8, K10, K11 closed).
- LIMITATIONS 24 (the traces: the nested reading fixes leg 2's rows; reading B rejected; A7 and D5 and their counts),
  28 (Fig. 4: E5), 29 (arrays' inter-copy rule; the two conflict pairs; the extension's tries), 30 (new:
  what `check` predicts, K13, K14), E (the two rc7-final items fixed; the newly refused inputs).
- REPRODUCTIONS, PRESET_PROVENANCE, DISSERTATION_NOTES § 3 item 7 (Fig. 4, the traces) and § 4 (the 8 µm start's
  sources; Fig. 4's open item), README (the tests block; exit code 6), EVALUATION § 3 (what each suite needs:
  `test_rc7` the photomasks, `test_figure_anchors` and `test_rc8` the run pickles), ALIGNMENT S1 (Xia's sentence on the
  taper), foundation/README (which glossary copy is v114), OPTIONS (the floor start's Q2 cell; the Fig. 4 rows).
- `tools/dissertation_figures.py` (new, outside the stamp): main-text figures from a finished build (the photomask
  with four zooms; the realised cutoff, the filtrate fraction and the model's scenarios; a provenance table) and the
  V13 figure.

## 1.0.0rc7 (final freeze 2026-10-19 18:39), code `52ba91eddb3f07b8`

rc7 was first frozen at code `4a624626f570918c` (2026-10-19 14:23; the section below). Three of A. H. Routt's decisions
after that freeze changed it. The foundation (App. A, App. B, glossary v114, spec v5) is not edited.

1. **The clearance guarantee.** 16:08: a layout figure showed a meander head drawn across another leg, and he asked
   "how can we ensure that this never happens?" His decision, "Fix in rc7 now", unfroze the code (register BUILD-2).
2. **Dinh's figure-aligned cutoff is fitted to the printed device length** ("Fit printed 75.8 mm"; register INV-1,
   revision).
3. **Gifford Fig. 4 keeps its 8 µm start, with a 20 µm folded-head variant shown beside it** ("Keep 8 µm + variant";
   register INV-2).

Every record was rebuilt at the final code and compared with the first freeze's
(`provenance/rc7_scripts/compare_frozen_final.py` → rc7's `logs/compare_frozen_final.{md,json}`, in the rc7 evaluation
archive): of the 53 records (51 builds
and the 2 arrays), 52 are the first freeze's: every number, and the photomask byte for byte where one is written (m14's
now under its `_FAILED` name; `wbc_xia_process_refused` is refused, as before); the 3 Xia presets whose input lost
`meander.rectilinear_rows` are among them. The one that changed is the Dinh figure-aligned refit (INV-1), as intended.
One record field changed with the guarantee on every build that runs the DRC: `validation.drc.file` names
`<name>_UNASSESSED.dxf`, the file the build audits before its verdict names it (LIMITATIONS section E; the comparison
reports these "identical apart from the DRC's file name"). The Fig. 4 variant `gifford2014_fig4_start51` was rebuilt
too, byte-identical. Against rc6 (rc7's `logs/compare_rc6_rc7.md`, regenerated; in the rc7 evaluation archive): 48 identical, 47 photomasks byte-identical, 1
changed with its input (Dinh depicted, square turns), 4 only in rc7, as at the first freeze.

All eleven suites pass at code `52ba91eddb3f07b8`, tests `dc3ca9329e04088f` (rc7's `logs/test_*.log`, in the rc7 evaluation archive): test_appendix_trace,
test_units, test_bend_fd, test_options_sheet, test_negative (63 checks), test_presets_match_replicas,
test_figure_anchors (42), test_lattice, test_trace_anchors (268; the Dinh figure-aligned mask retraced: only Fig. 2B's
transition panel moved at two decimals, its median 0.71 → 0.74 µm), test_rc7 and the new test_composition (96).

### Code (register BUILD-2)

| change | what it does | effect on a design that was clear |
|---|---|---|
| collision gate before writing (`build.collisions`, `build.plan_and_draw`, `build.LayoutCollision`) | The build draws the outline in memory and checks it before the photomask is written: every part (the inlet with the head, each leg, each U-turn, the outlet) at least the floor from every part it is not joined to; outlet ports separated; every outlet port at least 2 mm from the inlet port. An outline that fails is REFUSED (exit 5, tag `layout`). The record names the parts, the shortfall and the remedies (`pipeline.layout_remedy`; for a given paperclip spacing, the spacing the generator's own rule gives), with a figure, `<name>_refused_layout.png` (`renders.refused_layout`). No photomask is written. | none: the same checks, computed once and reused for the record |
| `opencif check` and `generate(check_only=True)` draw the same outline (`pipeline.layout_check`) | They refuse what a build refuses (rc7's first freeze planned without laying out, so it passed the spec that collided). `check_only='plan'` keeps the plan-only check, faster; `tests/test_lattice.py` uses it. `generate(check_only=…)` now plans the loaded spec, as a build does (an alias such as `process.wall_along_flow_um` was planned at its default, and `filtrate_to_retentate_ratio` raised a KeyError; code review 4, F9). | none |
| outlet-stretch extension (`spec.L_out_extra`: `layout._plan_paperclip`, the serpentine stretch, a straight layout's stretch, `layout.stretch_total` with a flare) | Where only the outlet meets the inlet end (a folded layout whose last leg ends beside its inlet with no loop there to clear, e.g. two nested legs, or a straight device too short for its ports), the stretch is lengthened by the ports' shortfall plus 200 µm (at least 500 µm a step, up to 12 times) and the design re-planned; its closure takes the longer stretch. It stops when a step does not help, when another pair of parts collides too, or when the longer stretch no longer closes the outlet; the gate then refuses the best outline drawn and the refusal lists the tries. Recorded as `design.outlet_stretch_extended_um` and a note. Not for a paperclip's ring outlet (the ring sets that stretch). With `legs: auto`, when every count collides, the fitter hands the build the first count whose outlet alone meets the inlet end. | none: applied only to an outline that collides |
| `pipeline.name_mask` | The build writes `<name>_UNASSESSED.dxf`; its verdict names it: a FAIL (exit 2) `<name>_FAILED.dxf`, kept for diagnosis and never under the deliverable name, any other verdict `<name>.dxf`; a re-assessment that clears a FAIL renames it back. Arrays likewise; `name_mask` does not replace an earlier `<name>.dxf`, which `pipeline.run` removes first but `compose.array` does not (the final document check, M1: LIMITATIONS section E; `tools/make_arrays.py` now clears the array's folder, and both shipped arrays were composed again). A composed array whose copies or manifold merge or disconnect, or whose product port joins a filtrate port, is not built, and no photomask is written. | m14 (a deliberate Q2 FAIL): its photomask is now `m14_…_FAILED.dxf`, byte-identical |
| `spec.conflicts_with_reproduction` | Refuses (spec error) a value the spec gives that a reproduction field replaces with a different one, which rc7's first freeze only listed in `defaults_used` (code review R5b): the spec keys mapped in `spec.SPEC_ENGINE` (e.g. the head's width or rows, a fillet), and (code review 4, F1) the arms of an arms head, the inlet's port gap, the outlet's stretch and ring, the tubing's resistances.  Not covered (the final document check, m5): `meander.enabled: false` beside `reproduction.head`, and `meander.narrow_um` beside `reproduction.w_arm` (LIMITATIONS 29). Also refuses a leg spacing given beside `reproduction.leg_Y`. | The three Xia presets that depicted the head carried `meander.rectilinear_rows: 9`, which the depicted head's 76 rows replaced; the value is removed and the head drawn is unchanged. rc6's archival variant specs `provenance/rc6_traces/variants/xia2016_rs7_depicted_{fig3,dfill}.yaml` and `rc5/xia2016_rs7_depicted_s181.yaml` carry the same value and no longer load; they are records of rc6's builds, which rc7's traces read and do not rebuild. |
| `layout.leg_linear_local` (a corner the sweep found) | Where a leg's last wall block runs past the leg's end into a U-turn, its upstream corners' arcs end on the band edges before the leg's end (an arc cut there left a 0.19 µm step against the turn's wall, which the corner-radius check FAILed: a short paperclip whose leg 1 holds only the head). | only a leg whose last opening ends closer to its end than the corner radius: none of the shipped devices (no leg of theirs ends on a wall-block row: checked on the 50 built records' U-turn rows and row kinds) |
| `tests/test_composition.py` (new suite) | 22 devices sampled from layouts, heads, inlets, outlets, tapers and start rules on four processes, at the plan stage, with 8 full builds and an array: every check answers; the generator's own placements are outlined and never refused for colliding parts, their outlet stretch lengthened only where the outlet meets the inlet end (two nested legs, +1803 µm); a leg spacing given too tight for a head (a reconstruction of the case above: the same process and head style, the legs 818 µm apart) is refused by the check, naming the collision and the remedies with the generator's spacing, and by the build, with no photomask written and the refusal figure; the 7 drawn builds pass every hard geometric check with the exit code and outlet stretch their check predicted (case 00 declares its Q2 FAIL); photomask names follow the verdict.  96 checks, ALL OK at the final code (`logs/test_composition_A.log`, `_B.log`: 30 and 32 min, one after the other).  The case itself (`provenance/rc7_scripts/clearance_original_case/`): `opencif check` and `opencif build` on its spec both refuse it at the final code (exit 5 in 50 and 66 s; 'inlet and head' and 'leg 2' overlap; no photomask; the refusal figure), where the first freeze drew and wrote it and FAILED it after 14 minutes.  Before code review 4's fixes the sweep found the corner above (case 12, a FAIL) and a 2-leg nested paperclip whose outlet met its inlet (the extension). | none |
| independent code review of the change (`CODE_REVIEW4.md`, in the reviews bundle; it reviewed code `0763a945e8f47c7b`) | Two major findings (contradictory-input pairs the conflict rule missed; the test's build condition), eight minor and two nits; all fixed before the freeze except the arrays' port spacing (F10), which LIMITATIONS item 29 states as not covered.  The fixes were not re-reviewed as code; the final document check (`FINAL_CHECK.md`, independent, read-only) checked them against the diff and found them as described, apart from its M1 (an array composed into a used folder keeps the earlier photomask under the deliverable name; LIMITATIONS section E) and m5 (two contradictory-input pairs the rule misses; LIMITATIONS 29). | none |

### Presets and variants

| item | change |
|---|---|
| `dinh2024_v1_figure_aligned` | D 4.29 → 3.67 µm, fitted so that the drawn device (the depicted head, the linear section, the inlet with its filter, the outlets) is the printed 75.8 mm long. Result: 75.71 mm; N 1093, 874 linear; exit 0. Dinh's 4.7 µm beads agree: on Gifford's measured soft edge their passage implies D 2.3–3.9 µm. Caveats: the figures' proportion is given up (the linear section starts at 41.5 % of the inlet → middle-outlet distance, against 46–58 % on Fig. 3A); 133 linear rows (21 % of the extraction) run above κ_g's calibration cap. Header comments, `paper.targets_fitted` and the description rewritten. Analyses: `provenance/rc7_scripts/postfreeze/`. |
| `xia2016_rs7_aligned`, `_depicted`, `_depicted_pullR` | `meander.rectilinear_rows: 9` removed (never drawn: the depicted head's own 76 rows are; the photomasks are byte-identical) |
| `provenance/rc7_traces/variants/gifford2014_fig4_fold20.yaml` (new) | Gifford Fig. 4 at a 20 µm floor, with spec v5's FOLD start (a sealed head), legs 1–2 at the generator's own spacing for the head. Built, traced and shown as a figure, not a reproduction (Gifford drew no head; the inset's gaps read 17.2 µm). Built at the final code: exit 0, N 3546 (a 59-row head), every check passing, no row's cutoff above the 8.3 µm beads; traced on the inset: the claim's median 1.10 µm, p90 2.51 µm (at its 20.5 µm gap, as the design preset's 1.10 / 2.51). |

### Documents

- LIMITATIONS item 29 (parts that collide: what the gate covers, the extension's limits, and what it does not cover: joined parts, ports apart from the 2 mm inlet–outlet rule, arrays' inter-copy port spacing, the sample), item 28 (Fig. 4's start: the beads' bracket and the fabrication sources), and section E (the corner the sweep found, fixed; the DRC's file name in the records, deferred).
- Tools (not in the code stamp): `tools/make_arrays.py` clears the array's earlier outputs before composing; `tools/package.py` ships `logs/compare_*` and `logs/rebuild_rc7.txt`; `tools/options_sheet.py`'s exit-3 text (the check outlines in memory).
- `provenance/rc7_scripts/clearance_original_case/` (new): the spec A. H. Routt's figure came from, refused at the final code by `opencif check` and `opencif build`, with the refusal record and figure.
- `provenance/rc7_scripts/compare_frozen_final.py` (new) and `compare_rc6_rc7.py`: a FAIL's photomask is read under its `_FAILED` name, and the DRC's file name is reported, not counted as a change.
- RESULTS and CONTACT_SHEET (regenerated from the final records), README, OPTIONS (regenerated), EVALUATION (eleven suites; the guarantee in the claims table), DISSERTATION_NOTES, REPRODUCTIONS, PRESET_PROVENANCE: the Dinh refit and the guarantee.
- The register: BUILD-2, the INV-1 revision and the INV-2 additions.
- The fleet's ±0.5 µm bias statistic on a headed device: up to 23 % of N (`dc_D3_F0p3`, `dc_D4_F0p3`, as at the
  first freeze), stated now apart from declared rows: the first freeze's 29 % and the refit's 24 % on Dinh
  figure-aligned are mostly its declared head rows, outside at nominal (the final document check, m4).  Under −0.5 µm,
  24 of the 33 headed designs (was 23: the refit now loses one seam row) lose one row's cutoff.

## 1.0.0rc7, first freeze (2026-10-19 14:23), code `4a624626f570918c`

From the critical re-check of rc6's four device packages and their register (`RC7_DISCREPANCIES.md`, delivered with the
packages), A. H. Routt's instruction of 2026-10-18 23:09 ("Regarding Dinh, rc7 should investigate it"), and the rc6
reviews' adjudication (its § 12-C items, here ADJ-C1–ADJ-C9 so as not to be read as App. B's claims C1–C9).  The foundation (App. A, App. B, glossary v114) is not edited.  Every
new drawing option is off by default (the passage report, ADJ-C9, is on, but it writes the record only), so every rc6
photomask is expected to rebuild byte-identical unless its preset opts in or the data it reads changed; the rc7
rebuild checks this against rc6 (§ "Rebuild against rc6", written from the rebuild's comparison).  With BUILD-1 the
masks match rc6's where rc6 built each spec in a fresh process, as `tools/run_batch.py` does.  The independent code
review of rc7 (codereview1, findings R1–R16) is answered in its own table below.

### Code

| item | change | default (rc6 unchanged) |
|---|---|---|
| BUILD-1 | `cifsched.clear_caches()`, called first by `pipeline.run`: the bend and depth-weight memo caches (keyed by rounded inputs) made a photomask depend on what the same process built before it (three coordinates of the Xia package's second photomask moved by 10⁻⁵ µm) | always |
| XIA-1 | `centre_channel.shear_formula`: the shear the constant-shear taper holds, `exact` (CIF_v2's C(α) Q/(a b²), spec v5 S16-b) or `parallel_plates` (6 Q/(a b²), so w_c ∝ Q_c): the form Xia's Fig. 3 panels d and e read, a departure from spec v5 S16-b (constant shear with C(α), S:456; "max γ̇ via C(α) [NEVER 6]", S:606), so reproduce mode only (Xia's text does not state the coefficient: the formulas it cites are Yang et al.'s, whose Eqn 3 the glossary identifies as the resistance series, G:2036–2038).  The record reports the true peak shear via C(α) (`validation.gates.shear_peak`, S:606) | `exact` |
| DINH-1 | `meander.turns: square` with `corner_radius_um`, `inner_corner_radius_um`: the head's quarter turns drawn as right-angle corners (Dinh Fig. 2B), the centreline through the corner point, plus a corner term c_corner w from 2-D solves (`provenance/rc7_scripts/corner/`: Hele-Shaw and Stokes, validated on rc6's sharp-L and exact-arc values, depth-weighted and banded as the junction terms; zero under `junction_band: none` and `bend_model: none`).  The head-capacity lengths (`Spec.L_taper`, `L_cap_w`, `R_of_L`) take the corner's 2 r.  The drawn-mask measurement takes the roundings' area off before converting a unit's area to length (the corner term counts it): found by the first square-turn build, whose head units all read 0.39 µm per corner long (Q1 0.33 % at row 222; after the fix 0.164 %, rc6's value).  The record carries the corner term (`model_scenarios.bend_model.square_corner`) | `arc` |
| GIF4-1 B | `layout.order: nested \| path` for the paperclip: legs side by side in path order (Y = 0, d, 2d, …, or a `reproduction.leg_Y` stepping one way), the loops alternating their turning direction | `nested` |
| unnumbered | `validation.gates.plateau_u_min_all_planned_rows`: the smallest planned row cutoff over its centre width, declared rows included (reported); a note when declared rows lie below the κ plateau while the device's D is on it | none |
| ADJ-C1 | `api.py` docstring: the two design targets are the only targets; depth, gap, posts, flow and fluid are model inputs | none |
| ADJ-C2 | `data/r6/rc6_values.json` note names the Xia arc value (values unchanged) | none |
| ADJ-C3 | Dinh presets: "row pitch" for the glossary's word | none |
| ADJ-C4 | `process.wall_along_flow_um` / `wall_cross_flow_um`: the glossary's names, aliases of `post_length_um` (g_p) / `post_width_um` (ℓ_g); both given must agree | none |
| ADJ-C5 | `fluid` section (`density_kg_m3`, `viscosity_Pa_s`, `rheology: newtonian`): the Reynolds and Stokes gates only (the schedule's resistance ratios do not depend on μ) | water |
| ADJ-C6 | `assessment.status_codes` (enumerated, `pipeline.STATUS_CODES`), `assessment.status.source_trace` (from `paper.source_trace`, a pointer: NAMED, or NOT NAMED without it; the build neither reads nor checks it, `tests/test_rc7.py` checks that each shipped pointer resolves to traces of its shipped mask), both printed by the CLI | none |
| ADJ-C7 | the record's `environment`: the test suite's and lock file's hashes, the installed versions against the lock, Python and platform (beside the code stamp, not in it) | none |
| ADJ-C8 | `defaults_used`: a value a reproduction block sets instead is named as the value drawn; derived values carry their effective value | none |
| ADJ-C9 | `validation.passage` (`checks.passage_report`): each measured channel class against D and 2 D, a warning note below 2 D; a report, not a check | on |
| DINH-7 | `meander.junction_radius_um`: the wall corners where each head unit's stem leaves the centre channel (the head gaps' entrance) rounded (`layout.junction_roundings`); drawing only, since the gap entrance coefficient c_ent is an external input | 0 |

### Code review of rc7 (codereview1, findings CR1–CR16; the register and the presets cite them as "code review R…")

| item | change |
|---|---|
| CR1 | square-turn corners below the head fillet are drawn, measured and modelled at the fillet (the head smoothing rounds every corner at it): `spec.to_engine` takes max(radius, r_fillet), notes a given radius below it, and lists the derived radii with their effective value |
| CR2 | `data/depicted.json['dinh']` re-derived for the square-turn head (`provenance/rc7_scripts/regen_depicted_dinh_square.py`: f +0.27 … +0.95 %, row 4 −1.31 %; the `xia` entry unchanged), and `depicted.py` derives it so; the depicted preset's head is now 219 rows (220 with the arc-derived fractions); the figure-aligned preset's head is written by `provenance/rc7_scripts/make_dinh_figure_aligned_head.py`; the depicted preset's length notes recomputed on this tree (`provenance/rc7_scripts/dinh_lengths/`: the schedule 0.11 mm shorter; F at the taper's end 0.5104) |
| CR3 | the passage report also compares every class with the largest particle a gap passes, min(D_realised_max, G), when rows are realised above D, with a warning (the Fig. 4 aligned start's 8 µm side channel against rows 1–3's 10.2 µm) |
| CR4 | the rc7 notes (passage, declared rows below the plateau) replace their earlier copies, so `reassess` does not repeat them |
| CR5 | `defaults_used` lists a value only where it applies (no meander values without a head, no inlet stem without a filter or an arm head), and names a given key that a reproduction block replaces (Xia: `rectilinear_rows: 9` given, 76 drawn) |
| CR6 | `plateau_u_min_all_planned_rows` over the rows extracting at least 5 % of the constant-cutoff fraction at their width (Xia's depicted rectilinear rows extract 2 × 10⁻⁷ of it), with `plateau_F_share_below_all_planned_rows` |
| CR7 | a square-turn unit's head capacity counts its corner term (`Spec.corner_term`); arc heads unchanged (LIMITATIONS 20) |
| CR8, CR9 | the figure-aligned Dinh preset no longer carries the depicted preset's length analysis; the Fig. 4 aligned header says rows 1–3 lie above D |
| CR10 | the traced presets name their rc7 trace records (`paper.source_trace`) |
| CR11 | `provenance/rc7_scripts/measure_p1.py`, `gifford_alignment.py` and `dinh_corners/` shipped with the tree (the presets cite them) |
| CR12 | behavioural tests in place of source checks where a fast one exists (the caches, `run` clearing them, the junction roundings), and new ones: radii below the fillet, every heading of a square corner, the passage report's idempotence, `reassess` twice, and rc6 invariance of arc-turn heads by fingerprints of the planner and the head geometry computed on rc6 (`provenance/rc7_scripts/fingerprint_head.py`) |
| CR13 | a square corner's outer cut reaches 10⁻³ µm beyond the walls, so no zero-width spike survives for the smoothing to grow (it refilled part of one rounding in 66 on Dinh's head) |
| CR14 | square-corner radii checked on the head's drawn unit widths: outer ≤ 0.307 w (the solved range), inner ≤ half the narrower of the unit and the wall; the record reports the corner table's deep and shallow values under `bend_model: none` and notes that the `junctions_ignored` scenario drops the corner term |
| CR15 | this header; README (the fluid, the rc7 options, `tests/test_rc7.py`); the record's absolute resistances (Pa s/m³) are labelled as water values |
| CR16 | check-only records carry the environment and clear the caches; the environment stamp closes its files and reads the lock as bytes; the status has a `reproduction` field in every mode; the status-code check survives `python -O`; `paper.source_trace` takes a string or a list of strings (anything else is noted and ignored, not refused); trace-code fixes; the `_plan_paperclip` docstring |

### Second code review and second foundation check (codereview2, foundation check 2)

| item | change |
|---|---|
| status code | `assessment.status_codes.source_trace` is `NAMED` / `NOT_NAMED` (was TRACED / NOT TRACED; the status sentence reads "NAMED: ..." or "NOT NAMED: ..."): the build cannot know the trace, which is made from its own photomask afterwards; `tests/test_rc7.py` checks that each shipped preset's pointer resolves to traces of its shipped mask (codereview2 N1) |
| `depicted.py` | its `__main__` compares with the shipped `data/depicted.json` and writes only with `--write`, the caches cleared before each entry and the `xia` entry kept from the file unless `--xia` (codereview2 N3: re-running it could otherwise move the Xia masks in the last bits) |
| peak shear | `validation.gates.shear_peak`: the peak wall shear via C(α) along the centre channel, over the taper, and, for the parallel-plate taper, its ratio to 6 Q/(a b²) (spec v5 S:606; foundation check 2) |
| small | the corner-radius refusals name `reproduction.r_fillet` when the fillet comes from it; the solved outer-rounding limit 0.307 is one constant (`cifsched.CORNER_ROUT_SOLVED`); the record's depth weight under `bend_model: none`; rows extracting f ≥ ½ (no cutoff) make the passage report compare with the gap itself (`realised_cutoff.rows_f_half_or_more`) |
| tests | `test_presets_match_replicas` builds the depicted Dinh replica with square turns (it failed on rc7: codereview2 N2); the BUILD-1 stub test builds into a temporary folder; rc7's vacuous "plan after another" check removed; the source-trace pointer test |
| `dinh2024_v1_figure_aligned` | D 4.5 → 4.29 µm: the figures' row count gives 4.28–4.79 µm at the 1-D F and 4.04–4.52 µm at the measured F, and the beads' soft edge favours the low end; 4.29 µm is the printed-length reading (642 rows), inside every range (foundation check 2) |

### Third code review (codereview3, `review7/codereview3/CODE_REVIEW3.md`)

| item | change |
|---|---|
| B1 | `gates.shear_peak.over_the_taper` counted the meander head's rows (drawn at w_c0) as taper rows, so every exact constant-shear device with a head recorded a false drift (Dinh as specified 5.9 %, m02 4.5 %, Xia 0.29 %; the true spread is ~3 × 10⁻¹⁵): now the linear rows above the planner's floor max(w_c,min, 0.3 w_c0); the parallel-plate Xia device's drift moved from 1.1747 to 1.1740 (+17.4 %), since the head's rows at w_c0 no longer count |
| B2 | `realised_cutoff.rows_f_half_or_more` counted every row without a cutoff, including sides at f = 0 (Xia: 7 / 8 / 6 rows; the true count 0): now the rows with f ≥ ½ on a side; no passage verdict changed (Xia's G was already the bound) |
| m1, m2 | the passage note names what set the largest passing size (the f ≥ ½ rows, or min(D_realised_max, G)), and no longer fails with `D_max` absent; `passage.rows_f_half_or_more` recorded |
| m3, n2 | `depicted.py`: `--write` refuses an entry that differs from the shipped file without `--force`; without `--xia` it says the `xia` entry is kept, not "equals" |
| n4 | `shear_peak.formula` says it is evaluated on the planned centre widths at the operating flow |
| M1 (tests) | the source-trace pointer test expands brace groups; the Xia presets' pointers are now lists of exact paths, the Fig. 3a whole-device records included (checked by their `dxf_sha256`) |
| not changed | n1: `shear_peak.taper` falls back to `'const'` on a namespace spec without `taper` (test namespaces only; every shipped spec carries it) |

### Presets

| preset | what | codes |
|---|---|---|
| `gifford2014_fig4_aligned` (new) | reading B: the printed design and Gifford's start rule at his 10 : 1 limit (8 µm), the four legs in path order (`layout.order: path`, 818 µm), leg 4's array 8.2 mm short (fitted to the neighbour's width), G 17.2 / g_p 28.6 µm (measured).  The rc6 Fig. 4 package's device (a 51 µm start) is rejected by the inset's beads (INV-2) and kept as `provenance/rc7_traces/variants/gifford2014_fig4_start51.yaml` | GIF4-1..3; INV-2 |
| `gifford2014_fig6_aligned` (new) | the Fig. 6 package's device: printed design with the text's 144 µm, 40 µm start (fitted to leg 1), G 17.1 / g_p 34.3 µm (measured) | GIF6-1..3 |
| `xia2016_rs7_aligned` (new) | rc6's fig3 variant with `shear_formula: parallel_plates`, D and F refitted to Fig. 3a's legs and panel g (D 6.4159 µm, F 0.94072) | XIA-1 |
| `dinh2024_v1_depicted` | `meander.turns: square`, outer 7 µm, inner 10 µm (Fig. 2B) | DINH-1 |
| `dinh2024_v1_figure_aligned` (new) | a two-cutoff reading of Dinh's figures: the head exactly as the depicted preset plans it (`reproduction.f_head`, 219 rows: 18 measured on Fig. 2A's inset and derived for the square-turn head, then 201 at the printed f*; `NM_force` 219; written by `provenance/rc7_scripts/make_dinh_figure_aligned_head.py`), the linear section at D 4.29 µm (fitted to Fig. 2A at the printed 75.8 mm length: 642 linear rows; through the closure at the printed seam and outlet widths and App. A R13(a)'s invariant lamina the figures' 514–642 rows give 4.28–4.79 µm at the 1-D F, 4.04–4.52 µm at the measured F; Dinh's bead data, App. B P-C7.2, put the cutoff below 4.7 µm, soft edge); the record (code 4a624626f570918c): N 860 (219 in the head, 641 linear), F 0.55964, exit 0, footprint 65.2 × 3.16 mm (array 53.1 mm, head 25.6 mm; the printed reading's 148.2 mm; Dinh's chip 75.8 mm, Fig. 2A at the lateral scale 60.6 mm); its two caveats, the scale readings behind its head and its D, and the linear rows beyond κ_g's calibration, are in DISSERTATION_NOTES § 3 item 7 and PRESET_PROVENANCE | INV-1 |

### Rebuild against rc6

Every spec (51) rebuilt at code `4a624626f570918c`, each in a fresh process (`review7`'s resumable batch over
`tools/run_batch.py`'s list, two lanes), and the two parallel arrays; compared leaf by leaf with rc6's records
(`provenance/rc7_scripts/compare_rc6_rc7.py`; its output in rc7's `logs/compare_rc6_rc7.md` and `.json`, in the rc7 evaluation archive):

- **46 specs built in both trees with an unchanged input, and the two parallel arrays (`m07_…_x3`, `m13_…_x4`): every
  value of every record equal to rc6's** (apart from the timings, stamps, titles, logs, notes and the `defaults_used`
  text, which the comparison skips, and the arrays' element code stamp), and every photomask byte-identical
  (45 DXFs and the 2 arrays; `wbc_xia_process_refused` is refused, as in rc6, and draws none).  The only
  record fields that differ are those rc7 adds: the status codes and source trace (ADJ-C6), the environment (ADJ-C7),
  the passage report (ADJ-C9), the fluid's rheology (ADJ-C5), the square-corner record, the plateau over all planned
  rows (CR6), the resistance note (CR15), the peak shear (S:606, foundation check 2), the f ≥ ½ row count
  (codereview2/3), and the engine spec's new fields at their rc6-equivalent defaults (`leg_order`, `mu_f`, `rho_f`,
  `r_corner_in`, `r_corner_out`, `r_junction`, `shear_formula`, `turns`).
- **`dinh2024_v1_depicted` changed with its input** (square head turns, DINH-1): N 2702 → 2701 (the head 220 → 219
  rows), F as drawn 0.559283 → 0.559226, 148.3 → 148.2 mm; its DXF differs.
- **Four presets only in rc7**: `dinh2024_v1_figure_aligned`, `gifford2014_fig4_aligned`, `gifford2014_fig6_aligned`,
  `xia2016_rs7_aligned` (REPRODUCTIONS, rc7 table), and the rejected Fig. 4 variant
  (`provenance/rc7_traces/variants/runs/gifford2014_fig4_start51/`).
- **Nothing changed without an input change.**
- Every build record (51; the two array records carry no environment) carries the code stamp `4a624626f570918c`, as
  the array records do, and the same test-suite hash (`environment.tests_sha256`
  `7929671ae732e849`; codereview2's N5, a hash that varied within the first rc7 rebuild, does not recur: the tests were
  not edited during this one).

The suites at this stamp: `test_appendix_trace`, `test_bend_fd`, `test_presets_match_replicas`, `test_options_sheet`,
`test_rc7`, `test_negative`, `test_units`, `test_lattice`, `test_figure_anchors` and `test_trace_anchors` (268 checks,
the rasters given) ALL OK (`docs/EVALUATION.md`); each suite's output is in `logs/`.

### Convergence review 1 (`review7/convergence1/CONVERGENCE1.md`): documents only

No code, preset value or record changed.  The fixes:

- **Disclosures added** (DISSERTATION_NOTES § 3 item 7 and § 4, LIMITATIONS 24, PRESET_PROVENANCE, REPRODUCTIONS, the
  presets' comments, the register): the figure-aligned Dinh preset's two caveats (its head and its D rest on different
  readings of Fig. 2A's scale; 307 of its 638 linear rows, 49 % of the section's extraction, above κ_g's calibration cap);
  partial tests and consistency checks labelled as such (Fig. 6's start, reading B's walls, XIA-1's panels); Xia panel
  f on the shipped aligned preset (INV-4: 39.9 / 36.0 µm); the rc7 traces' limits set with the results in view and their
  21 null-ratio flags; Fig. 6 aligned's declared entrance (rows 1–40, 9.9 D, 14.5 % of F); which rc7 presets are inside
  the validity envelope (Gifford Fig. 4 aligned only).
- **Numbers corrected**: the parallel-plate peak-shear drift 17.4 % (it was 17.5 % in seven places); the square-turn
  head's band, 121 rows = 111 head rows and the last 10; the fold pitch 2p = 94 µm at the 19 µm floor (110 µm is
  Dinh's drawn pitch with 27 µm walls); the fleet statistics recomputed on the 50 rc7 drawn designs (any-channel Re
  bound ≤ 10 on 38; outlet ±1 % moves 12–59 rows; the stretch-end band 6–58 rows on the 24 designs whose outlet changes
  width; the ±0.5 µm bias 17–223 and 82–404 gaps, up to 29 % of N; d(D_N/D)/d ln R_s,out −26 to −175).  Two numbers that
  were stale before rc7: the noise draw moves 0–2005 gaps (the docs said 0–2016; rc6's records give 2005 as well), and
  the outlet band's "up to 1.10 D" and "12 rows" (the records give 0.76 D and 10 rows).
- **The passage list** (LIMITATIONS 26): rc6's `passage_scan.py` read each design's planned closure widths, so it listed
  `xia2016_rs7_specified` at 11.36 µm, the closure width of a rectilinear transition row the mask draws at the
  transition's own width; rc7's report on the drawn mask reads 16.02 µm there, and adds `gifford2014_fig4_aligned`.
- **ALIGNMENT**: AG7b (rc7 supplies L_eq as the computed c_corner), AB25 (the `fluid` section) and a new § 6 for spec v5
  S16-b / S:606 and § 14.6.  The trace test's disposition text for L_eq still names only the exact bend solution: a
  stated, not a tested, disposition; editing it would change the tests hash every record carries.
- **The rc6 adjudication's § 12-A edits** (2026-10-18), planned for rc7 and not applied until now: A1–A4, A4b, A6, A8,
  A10–A12, A15–A17 (CHANGELOG, LIMITATIONS 24, EVALUATION's opening, PRESET_PROVENANCE, REPRODUCTIONS, OPTIONS through
  `tools/options_sheet.py`, the rc6 traces README and verification); A7 in part; A13 and A14 not done; B1–B3 dropped
  (register § 5 gives each status).
- **Release files** (rc7's; in the rc7 evaluation archive): `logs/batch_summary.txt` rewritten from the final legs (one line per spec), `logs/rebuild_rc7.txt`
  annotated, the suite logs copied into `logs/`, the tree-root `EVALUATION.md` refreshed from `docs/`,
  `provenance/rc7_traces/README.md` regenerated (the found-sample basis of the medians; Xia f's widened scale range).

## 1.0.0rc6 (2026-10-17), code `4bbdaabef1225a24`

From the reviews of rc5's micrograph traces (Reviewers A and B), the verification of those reviews
(`provenance/RC6_REVIEWS_VERIFICATION_2026-10-17.md` §§ 0–13) and the implementation and its independent checks (§ 14
there).  The foundation (App. A, App. B, glossary, frozen at v114) is not edited; its two numeric errata are recorded
for the next unfreeze (`foundation/ERRATA_NEXT_UNFREEZE.md`, D96 and D97).  Every rc6 change to the hashed package, the presets, the examples and the package data is
made by a script in `provenance/rc6_scripts/` (`apply_all.sh` runs them in order on the rc5 tree), so the changes can
be re-applied and reviewed as code; the tests, the tools and the documents are edited as text (the documents' trace
numbers taken from the records).

### Figure measurements: every scale bar read with its anti-aliased ends

rc5 read each bar at its full-intensity pixels only, which reads it short (Dinh's 50 µm bar 46 px for 46.906, 2.0 %;
Gifford's 300 µm bar 107 for 107.52, 0.5 %; Xia's Fig. 3a 5 mm bar 179 for 178.32).  rc6 re-reads every
figure-derived value, with two uniform rules: ℓ_g = the post columns' centre-to-centre separation − the printed centre
(free of the outline bias), g_p = Λ − G.  Where a bar cannot set a panel's scale (Xia's Fig. 3 panels b–g, whose bar is
panel i's; Gifford's Fig. 6, whose bar runs off the raster), the panel's printed widths set it and the bar is a check.
Scripts and outputs: `provenance/rc6_scripts/`; method, validators and traced figures: `provenance/rc6_traces/`.

| specs | change | rc5 → rc6 |
|---|---|---|
| Dinh V1 (both presets), `api.PROCESSES['dinh2024']`, `examples/dinh_targets_sealed_head.yaml`, `examples/platelet_serpentine_chip.yaml` | g_p (Λ 45.11 µm, six post-column readings on Fig. 2B panels 2–4, SD 0.04; ± 0.2 µm with the bar's reading) | 27 → 26.1 µm |
| Dinh V1 (both presets) | outlet flare onset (rc5's reading at the corrected bar) | 67 → 65.4 µm |
| Dinh V1 (both presets) | `paper.note`, `paper.length_bound`: recomputed at Λ 45.1 ± 0.2 µm by rc5's scripts ported to take the preset's own pitch (they reproduce rc5's outputs exactly at 46.0 µm); the source locator (Fig. 2 is on p. 917, the text on p. 916) | see the presets |
| Gifford Fig. 4 and Fig. 5 presets, the examples on Gifford's process, `api.PROCESSES['gifford2014']` | g_p; ℓ_g (Λ 45.82 µm on the Fig. 4 inset) | 25.6 → 25.3; 56.8 → 59.5 µm |
| Gifford Fig. 6 (both presets) | g_p; ℓ_g, measured on Fig. 6's own micrograph at its printed 125 µm centre between the post rows' tip lines (Λ 51.38 µm, three legs, synthetic calibration of the tip estimator; ≤ 52.01 µm by the cropped bar) | 32.2 → 30.9; 56.8 (Fig. 4's) → 70.8 µm |
| Xia RS-7 (all three presets), `api.PROCESSES['xia2016']`, `examples/wbc_xia_process_refused.yaml` | the scale of Fig. 3's panels b–g (their printed widths, post columns and loop curvature, not panel i's bar) | 1.953 → 1.798 µm/px |
| | g_p; ℓ_g; meander wall; head width | 40.25 → 35.8; 78.1 → 70.9; 41 → 35.8 µm; 2000 → 1840 µm |
| | Fig. 3a's layout (ports, ring, leg axes, offsets) | × 1.003827 |
| | the depicted head (`opencif/data/depicted.json`, regenerated) | 31 → 29 units, 65 → 76 rectilinear rows |
| Gifford Fig. 6 (both presets) | the layout, read on the micrograph after the author's review of 2026-10-17 (`apply_rc6_layout_presets.py`; PRESET_PROVENANCE, the convention): `reproduction.leg_Y` (legs 1, 3, 5, 7 on its rows; legs 2, 4, 6 from the turn walls' common centre, assuming semicircular turns) and `pc_last_short` | uniform 723 µm (turn radius 362 µm) → the micrograph's axes (the far-end turns concentric about a centre 1745.7 µm beyond leg 1, centre-line radii 1.75–2.92 mm); last leg 0 → 257 µm short |
| Gifford Fig. 4 (all three presets) | the leg spacing, read on the inset (the same review) | the generator's 969 → 818 µm |

New data files (hashed): `opencif/data/r6/rc6_values.json` (the values the engine reads), `opencif/data/r6/
figure_anchors.json` (the anchors `tests/test_figure_anchors.py` checks; the r5 file is kept as rc5's record).

### Code

* **E5 (Reviewer A): a drawing defect is a named FAIL, not a crash.**  `engine/measure2.py`: every measured width that
  enters r(w) (head units, gaps, side segments, the centre at the post tip and in the head, loops, the outlet stretch,
  flare and leads) goes through one guard; a section that finds no fluid where the plan puts the feature is recorded
  once per (row, feature, side) in `validation.features_found`, the design width stands in so that the measurement
  finishes, and `pipeline.assess` fails the build on the new hard check "every planned feature found on the drawing".
  A good drawing is measured exactly as before (the guard returns a positive width unchanged).
* **E5's cause: the fillet default.**  The head and outlet smoothing (a closing then an opening of radius
  `process.fillet_um`) erases channels and fills walls no wider than twice the radius; the 5 µm default did so at a
  10 µm floor.  `spec.py`: the default is min(5, 0.4 × the narrowest drawn feature, 0.3 × the narrowest head unit),
  the last keeping r_f / w inside the head T-junction term's fitted range; it is 5 µm for the 38 shipped specs that use
  the default (the other 9 set it explicitly, each below half their floor).  The narrowest drawn feature is the floor,
  the gap G, the post length g_p and, with a head, the meander wall, channel, arm and narrow widths.  An explicit fillet
  at or above half of it is refused before drawing (`SpecError`, naming the field the fillet came from and the feature);
  one above 0.3 of the narrowest head unit is noted in the record (the T-junction term is held at its fitted 0.3).
* `outlet.flare_start_um` default 67 → 65.4 µm (`spec.py`, `cifsched.Spec`); no shipped spec uses the default.  Records
  without a flare carry the unused field (`design.spec.flare_start`), which therefore reads 65.4.
* **Deferred from the rc5 reviews** (`provenance/RC5_REVIEWS_VERIFICATION_2026-10-12.md` § 8.2; applied by
  `provenance/rc6_scripts/apply_rc5_deferred.py`): (1) the outlet stretch-end term's domain, which every record states
  (`design.outlet.stretch_end.table_domain`); a design beyond a solved 2-D table's end (steps R 8, flares R 4) is
  refused in strict mode (`layout.OutletOutsideTables`, exit 5, the refusal naming the outlet: at the plan stage, or,
  for a ring, whose lead is re-solved while drawing, before the mask is written) and built with the end values held in
  evaluate and reproduce modes, the verdict of a passing or reproduced build saying so (rc5 clamped silently); no
  shipped design comes near the tables' ends (the largest width ratio is 3.87, `wbc_paperclip_teardrop`); (2) `reproduction.L_out` with
  `outlet.flare_deg` is a spec error (rc5 ignored it on a linear layout); (3) the fail-safe refusal's message names the
  last side channel's width.  Not taken: (4), the lead arc's interaction in c for teardrop leads (optional in the
  report), and (5), a DXF-based mask check in the archive (a test option).
* `engine/devspecs.py`: the first-generation factories `dinh()`, `xia()`, `gifford_bead()`, `gifford_plt()` removed,
  because they carried CIF_v2's reference-process posts and nothing used them; the replica factories read rc6_values.json.
* Comments and docstrings brought up to date (`cifsched`'s transition row length and ring defaults, `devspecs.xia_fig`'s
  docstring, `pipeline.outlet_precheck`'s).  `__version__` 1.0.0rc6.
* From the independent review of this code (`provenance/RC6_REVIEWS_VERIFICATION_2026-10-17.md` § 14.4;
  `provenance/rc6_scripts/apply_rc6_review_fixes.py`): `devspecs` imports under `python -OO` (rc6's docstring line
  raised on a stripped docstring); `validation.features_found` counts the missing features by kind (`by_kind`; the list
  keeps the first 20); the stretch-end refusal names the outlet style, as every other outlet refusal does.

### Tests

* `tests/test_figure_anchors.py`: reads `data/r6/figure_anchors.json`; Gifford pitch tolerance ± 0.2 µm; ℓ_g on the
  Fig. 4 presets; Fig. 6's pitch and ℓ_g against `provenance/rc6_scripts/gifford_fig6_measurement.json` (and, with
  `OPENCIF_PAPER_RASTERS`, the raster's hash against the measurement's); the Dinh length notes against
  `provenance/rc6_scripts/` at the presets' own pitch (23 checks need no build, from 17).
* `tests/test_negative.py`: 12b: `reproduction.L_out` with a flare refused; the stretch-end tables cut at R 1.5 so
  that a shipped design (m15) reaches them: strict refuses by name before drawing, evaluate plans and its accepted plan
  carries the status (the ratio and the table's end).  13 (E5, on
  Reviewer A's fresh device u05: 12 µm gaps, a 10 µm floor, 60 µm deep, the narrow head): an explicit 5 µm fillet
  refused by name, the default 3 µm; rc5's 5 µm forced past the spec rule: drawn and measured to the end, the vanished
  head units and gaps listed by name and the build FAILs on that check (rc5 crashed); at the default the same device
  is drawn with every feature found and its topology as planned; one gap filled on m15's drawing named by the
  measurement, the clean drawing as control.  Tests 8 and 10 print their pitch-dependent rates from the spec.  Test 5's
  plateau seed (a 3 µm floor on `examples/gifford_physics_linear.yaml`) lowers that example's explicit 2 µm fillet to
  1 µm with the floor: rc6's fillet rule refuses 2 µm at a 3 µm floor as a spec error before the plateau gate is reached
  (found when the full suite first ran on rc6).
* `tests/test_presets_match_replicas.py`: `exp_domain_refuse` (the strict-mode switch `spec.to_engine` sets from
  `gates.mode`) is not device geometry and is ignored with the other operating fields.
* `tests/test_units.py`: 12b checks the flared stretch against the spec's own default `flare_start` (rc5: a literal
  67); 13: the stretch-end domain at and beyond each table's end, and its refusal only when strict.
* `tests/test_trace_anchors.py` (the micrograph traces), rewritten for rc6: five layers (the code stamp; the inputs,
  with Gifford's Fig. 6 measured again from its raster; the shipped masks' panels (Dinh p1–p4, Gifford Fig. 4 and
  Fig. 6, Xia c–g and Fig. 3a), each re-registered, reproducing its anchors and judged against its limits; the variants'
  declared fields; the variants' traces); a panel or input missing from the anchors is a named missing input.

### Micrograph traces and whole-device figures (`provenance/rc6_traces/`; nothing there is hashed)

Every micrograph panel of the three papers that shows post rows, loops or turns, and Dinh's head (Fig. 2B, first panel), is traced with the shipped rc6 mask at the
panel's own scale, registered twice (features and chamfer) and judged against its nulls and a scale scan; Xia's head (Fig. 3b) is measured (the depicted head's unit pitch, width and reach are read on it) but not registered; Xia's Fig. 3h–i are not traced; Dinh's Fig. 3C is left out (its bar disagrees with Fig. 2B's: `provenance/RC6_REVIEWS_VERIFICATION_2026-10-17.md` § 3.4); the method,
validators and results are in `provenance/rc6_traces/README.md`, the figures (600 dpi, captioned; residual maps) in
`figures/panels/` there, every record in `results/`.

* Dinh Fig. 2B, all four panels: the mask's pitch is the image's (V5 along 0.996–1.000 on the linear and outlet panels);
  medians 0.44–1.29 µm; the outlet's widths match the print to ≤ 0.5 µm.
* Gifford Fig. 4 inset (leg 4): median 1.21 µm, V5 1.001 / 1.001; the three Fig. 4 presets take the inset's leg
  spacing (818 µm; the generator's rule gave 969), so the neighbouring leg's posts lie on the image's (consistency: the
  spacing is read on that inset).
* **Gifford Fig. 6** (new): the micrograph's four legs identified as 7, 5, 3, 1 (flow direction and the side widths at
  their ends, read without the mask); each registers on its own map (medians 1.4–1.9 µm), its along-flow pairing checked
  row by row (END_CHECK: the last pill and the mask's last post each one gap before the divider wall).  After the
  author's review (2026-10-17) the two Fig. 6 presets take the micrograph's layout (`reproduction.leg_Y`,
  `pc_last_short`: legs 1, 3, 5, 7 on its rows, 588 µm apart; legs 2, 4, 6 from the 16 turn walls' common centre,
  `code/gifford6_turns.py`, assuming semicircular turns): one map registers the four legs and their turns (the turn walls
  a median 1.4 µm from the mask's edges, against 23.6 µm with the generator's uniform layout, whose turns were too
  tight -- consistency, the axes being read on those walls; the turns' start along the flow, which they do not set, lies
  5.7 µm downstream of the image's).
* Xia Fig. 3 panels c, d, e, f and g at 1.798 µm/px (panel f at its own pitch); panel e recorded, not claimed (leg 3
  narrower than the constant-shear taper puts it).
* **Xia Fig. 3a, the whole device** (`xia_fig3a_device.py`): every part registered after the device is posed at its top
  (V11) and every traced panel against its box (V12): the shipped mask's far loops and outlet ring lie 5.69 / 5.41 mm
  short of the drawing's (its legs hold 1079 rows where Fig. 3a's hold 1190); the variant with D 6.5315 µm
  (`variants/xia2016_rs7_depicted_dfill.yaml`, the cutoff the drawing implies; F unchanged) registers end to end: far loops 0.09 mm, outlet ring 0.19 mm from the drawing's.  And the variant with F and D both read off Fig. 3 (`variants/xia2016_rs7_depicted_fig3.yaml`: F 0.94072 from panel g's last side channel, D 6.5705 µm from Fig. 3a's legs; `code/xia_fig3_fit.py`) matches panel g's wall (V5 across 1.005) and registers end to end (far loops +0.09 mm, outlet ring -0.19 mm); panel f, a quarter of the way down leg 7 (recorded at its own scale, not claimed), then reads its side channels narrower than drawn: if its scale holds, the device's side channels widen faster along the last leg than the model's (LIMITATIONS 24 d).
* Whole-device figures (`full_device/`): every polygon of the shipped DXFs with the micrographs placed by their stored
  maps (V9 round trip, V10 one scale), Xia's Fig. 3a under the device; the four shipped devices at one scale in
  `full_device_composite.png` (4 µm/px, assembled from their strips: `code/composite_from_strips.py`);
  `code/run_full_device.py` regenerates the whole folder for any tree.
* For a third-party reviewer: `provenance/rc6_traces/EVIDENCE.md` (`code/evidence_index.py`, from the records): figure
  by figure, the preset, the zoomed trace, the whole-device figure, the numbers, and what is tested against what
  agrees by construction.  `provenance/rc6_scripts/passage_scan.py`: every record's narrowest straight side channel
  against its D (LIMITATIONS 26).
* `tools/package.py`: a zip is closed before a file would take it over the 27 MB limit (rc5 checked after each file, so
  one large figure could overshoot); the trace figures and whole-device figures go to their own `opencif_traces_N.zip`;
  the suites' logs ship as `logs/test_*.log`.
* Fixed in the traces' own code (not hashed): `xia_arc_scale.py` used rc5's loop centres and failed on rc6's mask (it now takes them from the record's leg axes; panel d's arcs paired with rc6's walls give 1.792 ± 0.004 µm/px, re-anchored; the image radii, and so the rc6 scale's curvature term, are unchanged); `full_device.py` ran out of memory on Xia's devices (the underlay is cropped to each view, each figure is freed, the PDF raster resolution is a setting; each device is drawn in its own process).
* Found and documented: the side channels at Gifford Fig. 6's leg ends read 5.4–10.0 µm per wall wider than drawn
  (per leg 5.5 / 5.6 / 6.6 / 9.6 µm, growing as the channel narrows); a paperclip's
  last leg takes the leftover rows, so `reproduction.pc_last_short` moves its end in steps of about n − 1 rows
  (LIMITATIONS 24, 25).

### Records

All 47 shipped specs rebuilt at code `4bbdaabef1225a24` (`tools/run_batch.py`, logs in `logs/`), the two arrays composed from their rebuilt elements (`tools/make_arrays.py`: `m07_serpentine_const_pack_filter_separate_x3`, `m13_paperclip_const_nomeander_sealed_fine_x4`) and the 5 CFD-kit windows cut again (`tools/cfd_kit.py`).  Exit codes: 35 × 0, 1 × 2, 10 × 3, 1 × 5 (rc5: 35 × 0, 1 × 2, 10 × 3, 1 × 5); no exit code moved.  Against rc5 (`provenance/rc6_scripts/compare_rc5_rc6.py`, every leaf of every record): 30 identical apart from timings and stamps, 19 changed with a change of their input (the inputs rc6 re-measured on the papers' figures), 0 changed without one.
