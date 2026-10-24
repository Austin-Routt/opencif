# Provenance: the scripts behind values the package takes as inputs

Nothing here is imported by `opencif/` (it is outside the code hash).  These are the scripts that produced numbers the
package hard-codes or the docs quote, so an evaluator can see and re-run how each was obtained.

What 1.0.0 distributes: `junction_fits/`, `figure_measurements/`, `gifford_fig5_d50.py`, `rc4_scripts/`,
`rc5_scripts/`, `rc6_scripts/`, `rc6_traces/`, `rc7_scripts/`, `rc7_traces/`, `rc8_scripts/`, and
`RC8_DISCREPANCIES.md`.  The sections below also describe the review reports, independent checks, diffs, report
figures and superseded records of the release candidates (`RC*_REVIEWS_VERIFICATION_*.md`, `REVIEWS_2026-10-09.md`,
`rc3_independent_check/` to `rc5_independent_check/`, `RC*_DIFF.patch`, `RC5_REVIEW_EDITS.patch`,
`rc5_report_figures/`, and `runs_f6c8/`), which are not distributed; they are kept in the author's archive.  The RC8
module `cifsched` is `opencif/physics/schedule.py` in 1.0.0, and `measure2` is `opencif/checks/measure.py`
(`docs/CHANGELOG.md`, 1.0.0 item 1).  In `rc4_scripts/` and `rc5_scripts/`, the default folders of the release
candidates' trees are written relative to the working folder.

## junction_fits/: the 2-D fits behind the junction terms (`opencif/physics/schedule.py`, the RC8 module `cifsched`)

Standalone (numpy, scipy).  Each `*_conv.py` refines the grid and prints the converged values; the `.log` files are the
runs that set the constants.  `r3_junction2d*` print each bend's junction part as a fraction of its developed
correction; since the code removes κ_j·w of arc length from that correction, κ_j = (that fraction) × (arc length r_m θ) / w
(e.g. the Dinh hairpin U-turn, deep: 0.16 × π·27.5/28 = 0.49); over the converged geometries this gives 0.46–0.52 deep
(the code uses 0.50) and 0.33–0.35 shallow (0.33).  `r4_tjunc2d_conv` prints c_T directly (extrapolated
deep −0.281, shallow −0.465); `r4b_fillet_conv` prints c_T with fillets, whose change per unit r_f/w is S_T.

| constant in cifsched.py | value | script | log |
|---|---|---|---|
| `KAPPA_J_DEEP`, `KAPPA_J_SHALLOW` (a bend's straight-arc junctions remove kappa_j w of the developed correction) | 0.50, 0.33 | `r3_junction2d.py`, `r3_junction2d_conv.py` (deep: streamfunction-vorticity planar Stokes; shallow: Hele-Shaw) | `junction2d.log`, `junction2d_conv.log` |
| `C_T_DEEP`, `C_T_SHALLOW` (a head T-junction is shorter than the centreline convention by c_T w) | -0.28, -0.47 | `r4_tjunc2d.py`, `r4_tjunc2d_conv.py` | `tjunc2d_conv.log` (deep, stem 2: extrapolated -0.281; shallow -0.465) |
| `S_T_DEEP`, `S_T_SHALLOW` (fillet term, per unit r_f/w) | -0.83, -0.33 | `r4b_fillet.py`, `r4b_fillet_conv.py` | `fillet.log`, `fillet_conv.log` |

These are fits to 2-D solves in two limits, placed between them by a depth weight computed in the package from the
exact curved-duct solution: the weakest physics in the pipeline (LIMITATIONS 8).  Every record's model band re-solves
the drawn device with the terms at both ends of their range.

## figure_measurements/: the presets' measured (code M) values

Scripts that measured dimensions on the papers' figures (Gifford 2014 Fig. 4 inset; Dinh 2024 Fig. 2; Xia 2016 Figs.
2a, 3a, 3g).  Their outputs are the JSON files the package ships in `opencif/data/r4/` (byte-identical).  The scripts
read figure images extracted from the papers' PDFs, which are not redistributed here; each script's docstring names the
figure, the crop and the scale bar it used.  `docs/PRESET_PROVENANCE.md` codes every preset value.

## gifford_fig5_d50.py (+ .out): what the long Fig. 5 devices retain at their design cutoff

This project's reading of Gifford 2014 Fig. 5: the retained fraction R = E (1 − F) of each bead size, R(D) at each
device's design cutoff D, and the 50 %-retention size d50, under three interpolants and three values of the
full-retention ceiling (the duct law's 23.5×, the caption's 24.4×, the text's ~23 : 1).  At D the four devices retain
R(D) = 0.39–0.42 / 0.44–0.54 / 0.38–0.48 / 0.17–0.30 of the beads; d50/D = 1.09–1.14 / 0.98–1.03 / 1.01–1.04 / ≥ 1.10.
κ_g was fitted on an onset-of-loss threshold ("consistently maintained") on several-centimetre test arrays, so these
numbers describe cumulative loss in a long device; they are not a test of κ_g at its ±3.4 % level.  The duct-law
ceiling gives R > 1 for 5A and 5B's 10 µm beads, so the data favour the caption's ceiling.  (The "earlier, unshipped
analysis" figures 1.135 / 1.006 / 1.027 / ~1.108 are reproduced to 0.001 with linear interpolation by the caption's
ceiling for 5A–5C and the text's for 5D.)  New evidence, not part of App. A, App. B or the glossary.  Run from the repository root.

## RC1_REVIEW_VERIFICATION_2026-10-10.md

The author's verification, claim by claim, of the 2026-10-10 review of rc1 (verbatim quotes with locators, the
computations, what the review missed and got wrong, the rc2 edits and the rc1 → rc2 record comparison).

## RC2_REVIEWS_VERIFICATION_2026-10-10.md and rc3_independent_check/

The author's verification, claim by claim, of the two reviews of rc2 (Reviewer A; Reviewer B's alignment review and
ledgers): verbatim quotes with line locators, the computations, what each review got right, missed and got wrong, the
rc3 changes and the rc2 → rc3 record comparison.  `rc3_independent_check/` holds the three passes of an independent
checker over that report (33, 9 and 6 findings, all fixed in it).

## RC3_REVIEWS_VERIFICATION_2026-10-10.md, rc4_independent_check/ and rc4_scripts/

The author's verification, claim by claim, of the two reviews of rc3 (Reviewer A: its review of the rc3 verification and
its final review with a 37-case generator sweep; Reviewer B: an independent audit with claim and quote ledgers): verbatim
quotes with locators, the computations and measurements, what each review got right, missed and got wrong, the rc4
changes and the rc3 → rc4 record and mask comparison.  `rc4_independent_check/` holds the independent checker's passes
over that report.  `rc4_scripts/` holds the scripts behind its numbers; each takes the tree (or trees) on its command
line except `build_lattice.py`, whose paths are set at its top, and `fill_blat.py`, which reads `blat/` beside it:
`run_case.py` (one generator case, as the reviewer's sweep), `diag_ring2.py` (the ring outlet's lead target on rc3),
`step_sens.py` (the last rows' sensitivity to an unmodelled side-outlet resistance, on a preset's design network with its
drawn outlet), `build_lattice.py` (full builds of the 21 plan-stage lattice cases that plan inside the envelope and one Xia-process
case outside it, with its log and per-case results in `blat/`; `fill_blat.py` turns them into the report's §6.3 table), `ob_summary.py` (the outlet-balance numbers over a tree's records; output `ob_summary.out`),
`compare_rc4.py` and `deepdiff_tally.py` (the rc3 → rc4 comparison, given the two `runs/` directories: the first by
named fields, verdicts and notes, the second every changed leaf tallied by path; outputs `compare_rc4.out`,
`deepdiff_tally.out`).  `RC4_DIFF.patch` is the full rc3 → rc4 diff of the code, specs, tests, tools, docs and
provenance scripts (`diff -ruN`; records, logs, PDFs, generated tables and this round's outputs excluded).

## RC4_REVIEWS_VERIFICATION_2026-10-11.md, rc5_independent_check/ and rc5_scripts/

The author's verification, claim by claim, of the two reviews of rc4 (Reviewer 1: its verification of the rc4 report
and its implementation of the flared lead with 2-D Stokes solves; Reviewer P: a primary-figure audit with a claim
ledger, a quotation ledger and a patch plan): verbatim quotes with locators, the computations and measurements, what
each review got right, missed and got wrong, the corrections to the rc4 report, the rc5 changes and the rc4 → rc5
record and mask comparison; its two figures are in `rc5_report_figures/`.  `rc5_independent_check/` holds the
independent check of that report, four passes, each by a fresh checker that had not seen the work (pass 1's brief and
the four passes' findings, with the small scripts passes 2–4 wrote; their paths point into the author's scratch area and
are kept as written).  That report's statement that Dinh's length conflict holds at the presets' F only (none at the
measured F) is corrected by `RC5_REVIEWS_VERIFICATION_2026-10-12.md` (the every-row bound counts gap rows only; a
progressive taper conflicts at the measured F as well, at the measured pitch).  `rc5_scripts/` holds the scripts and outputs behind rc5's new numbers; each takes its inputs on its command line:

| script | what | output |
|---|---|---|
| `expfem.py` (+ `table.py`, `table2.py`, `conv.py`) | 2-D excess of the outlet stretch's end, Taylor–Hood P2/P1 finite elements (scikit-fem) on graded tensor meshes: the deep limit (planar Stokes, by the dissipation integral) and the shallow limit (Hele-Shaw, Laplace for p); one-sided steps, centred steps and one-sided linear flares; `conv.py` refines the mesh | `step_table.out`, `sym.out`, `sym_24_8.out`, `one68.out`, `flare_tables.out`, `tab_*.json`, `coefficients_all.json` (the values in `cifsched.STEP_ONE_DEEP`, `STEP_SYM_DEEP`, `FLARE_DEEP`, `FLARE_SHALLOW`), `conv.out`, `step_24.out` (the one-sided R 2.4 step on a finer mesh), `control.out` (the straight-channel control) |
| `check_tables.py` | the tables in `cifsched` against the solver outputs above (and the shallow steps against the closed form) | `check_tables.out` |
| `measure_dinh_flare.py` | Dinh 2024 Fig. 2B fourth panel: the outer walls' straight direction, flare angle and onset, on the publisher-embedded raster (`pdfimages -png -f 5 -l 5` of the Lab Chip PDF, image 0, 982 × 1940 px; not redistributed), every pick printed | `dinh_flare_measurement.json` (shipped also in `opencif/data/r5/`) |
| `measure_dinh_pitch_gap.py` | the same panel: post pitch and the last gap from the pill outlines | `dinh_pitch_gap_measurement.json` (also in `opencif/data/r5/`) |
| `overlay_outlet.py` | the drawn outlet end of a run pickle on the panel, registered on the last gap's centre and the channel axis only | the report's overlay figure |
| `dinh_length_bound.py` | the least row count and length the 1-D model needs for Dinh V1, at three readings of F, counting every row or the linear rows after the reconstructed head | `dinh_length_bound_depicted.json`, `dinh_length_bound_specified.json` |
| `rc4_step_excess.py` | rc4's unmodelled stretch-end excess on every rc4 record (rc5's coefficient on rc4's widths), as a share of R_s,out with the record's own derivative | `rc4_step_excess.json`, `.out` |
| `rc4_step_effect.py` | the same excess added to rc4's Dinh V1 network (the design network with the drawn outlet resistances), re-solved: rows leaving D ± 1 %, rows losing their cutoff, F | `rc4_step_effect_dinh.json` |
| `compare_rc5.py`, `deepdiff_tally.py` | the rc4 → rc5 comparison of the two `runs/` directories | `compare_rc5.out`, `deepdiff_tally.out` |
| `flare_rules_check.py` | the flared outlet's stretch rules on layouts no shipped spec builds: three examples with a flare added, planned (not drawn) on a separate outlet, and on a serpentine and a paperclip at one and three legs | `flare_rules_check.out` |
| `fig_coeffs.py` | the report's figure 1: the stretch-end coefficients from `cifsched`'s tables, with rc4's 21 changes of width from `rc4_step_excess.json` | the figure |
| `rebuild_diff.py` | the final rebuild against the first full rc5 rebuild (code `f6c8eba91ec1f316`, whose records are kept in `runs_f6c8/` as JSON with their DXFs' sha256 and `cfd_windows.sha256`): every changed record field, and every DXF (the CFD windows included); run as `python3 rebuild_diff.py runs_f6c8 ../../runs` | `rebuild_diff.out` |
| `band_rc5.py`, `stats_rc5.py` | the same two `runs/`: the model band per design (rows, the last rows among them, the largest change, rows losing their cutoff), and the record statistics the docs quote (outlet balance, aspect ratio, footprints, Re) | `band_rc5.out`, `stats_rc5.out` |

`RC5_DIFF.patch` is the full rc4 → rc5 diff of the code, specs, tests, tools, docs and provenance scripts (`diff -ruN`;
records, logs, PDFs, generated tables and this round's outputs excluded).

## RC5_REVIEWS_VERIFICATION_2026-10-12.md and the scripts added to rc5_scripts/ after the rc5 reviews

The author's verification, claim by claim, of the two reviews of rc5 (Reviewer 1: a verification of the rc5 report,
with an independent finite-difference re-solve of the stretch-end term; Reviewer P: an independent critical review
with a claim ledger and an edit plan): verbatim quotes with locators, the computations below, what each review got
right, missed and got wrong, the edits made (documentation, the Dinh presets' notes and their two records, a test's
exit status, provenance; the hashed package is unchanged) and those deferred to rc6.  New in `rc5_scripts/`, each
taking its inputs on its command line, run in the environment `requirements-provenance.txt` pins (Python 3.11.15;
scikit-fem 12.0.2, numpy 2.4.4, scipy 1.17.1, pillow 12.2.0, matplotlib 3.10.9):

| script | what | output |
|---|---|---|
| `regen_fem_tables.py` | every entry of `cifsched`'s stretch-end tables regenerated from scratch with `expfem.py` on the meshes they were made with, against the shipped values and the conformal closed form (one worker: two ran out of the 6 GB memory limit beside another solve) | `regen_fem_tables.out`, `.json` |
| `rounded_flare.py` | Dinh V1's flare as drawn (25 µm tangent arcs at both kinks) against the sharp-cornered flare of the tables, two meshes | `rounded_flare.out` |
| `step_bend.py` | a one-sided step followed at once by the teardrop lead's outward quarter arc (centreline radius 2 lead widths) against the same step with a straight lead, four geometries per width ratio, both limits, two meshes (`HMIN`, `HMAX`) | `step_bend.out`: three runs, W 2.4, 2.9 and 3.87 at hmin 0.005 / hmax 0.08; W 2.4 and 2.9 at 0.0025 / 0.06; W 3.87 at 0.0035 / 0.07 (at 0.0025 / 0.06 it ran out of the 6 GB memory limit) |
| `band_coverage.py` | on every record with a stretch-end term: the band's span against the deep–shallow range, and on the teardrop steps the c that `step_bend.out`'s interaction gives | `band_coverage.out`, `.json` |
| `dinh_schedule_at_F.py` | the Dinh presets' own schedules (head, assumed constant-shear taper) planned at F 0.5592, 0.5105, 0.505 and 0.472: rows and length | `dinh_schedule_at_F.out`, `.json` |
| `dinh_taper_alternatives.py` | the Dinh row count under `dinh_length_bound.py`'s premises (the lamina the same at every gap) for three taper shapes; the constant-equivalent extraction the rows need to fit the printed length; its sensitivity to the pitch's ±1 µm and to crediting the presets' heads; and, for comparison, readings that hold D itself constant with the frozen Hill curve below the κ plateau | `dinh_taper_alternatives.out`, `.json` |
| `anchors_scenarios/` | `tests/test_figure_anchors.py` on five simulated trees (no pickles; source anchors only; no provenance; one pickle missing; one preset value changed), with a README | the logs and `exit_codes.txt` |
| `rc6_patch1_trial/` | the verification report's rc6 patch (1) (a table-domain status with a strict refusal) applied to a scratch copy of the package and run with the table end forced lower, so the refusal path runs on shipped designs; `opencif/` is not changed (README inside) | `outputs/` (the CLI outputs and `status.txt`) |
| `measure_dinh_centre.py` | Dinh Fig. 2B fourth panel (the raster `measure_dinh_flare.py` reads, sha256 a21e2e3ce0710e32410ba4941f348c6a1ef561b748f569e24decce1dea9fbbe6): the centre channel's edges after the last gap, every pick, with a calibration between the post columns; the outer walls' straight-part tilt over six row ranges and the flare angle from each | `measure_dinh_centre.out`, `dinh_centre_measurement.json` |

`RC5_REVIEW_EDITS.patch` is the diff from rc5 as it was reviewed to this state (`diff -ruN` of the README, docs,
presets, tests, and the provenance scripts and READMEs; records, logs, PDFs, script outputs and the verification
report excluded).  The package itself
(`opencif/`), the tools and the examples are unchanged.

## RC6_REVIEWS_VERIFICATION_2026-10-17.md, rc6_scripts/ and rc6_traces/

The verification of the reviews of rc5's micrograph traces (Reviewers A and B; §§ 0–13, as delivered on 2026-10-16)
and rc6's implementation with its independent checks (§ 14).  rc6 changes the package only through scripts, so that
every change can be re-applied to the rc5 tree and reviewed as code: `rc6_scripts/apply_all.sh RASTERS` runs, in order,
the measurements, the code edits, the deferred rc5 items, the E5 guard, the replica factories, the review's fixes, the anchors, the depicted
heads, the presets, the Dinh length computations and notes, and Gifford Fig. 6's measurement and presets.  The rasters
are the publishers' figures, extracted with `pdfimages -png` from the project PDFs (sha256 in each output); the six that
the scripts read are shipped in `notebooks/figures/` (`NOTICE.md`).  Each script takes its inputs on its command line; Python 3.11 with `requirements-tested.txt`.

| script | what | output |
|---|---|---|
| `measure_rc6.py` | every figure-derived value re-read: the bars with their anti-aliased ends; Dinh's row pitch on Fig. 2B panels 2–4 and the flare onset at the corrected bar; Gifford's Fig. 4 inset (pitch, l_g by the centre-to-centre rule); Xia's panel scale (printed widths, column separations, loop-4 curvature), pitch, g_p, l_g, the head's unit pitch (wall) and envelope, and Fig. 3a's layout at its corrected bar | `rc6_measurements.json`; `opencif/data/r6/rc6_values.json` (hashed) |
| `measure_gifford_fig6.py` | Gifford Fig. 6's micrograph at its printed 125 µm centre: three legs, tip lines with a synthetic calibration, the cropped bar's bound, the end blocks as a cross-check | `gifford_fig6_measurement.json`, `out/measure_gifford_fig6.out` |
| `apply_rc6_code.py`, `apply_rc5_deferred.py`, `apply_rc6_e5_guard.py`, `apply_rc6_devspecs.py` | the hashed code changes (CHANGELOG): version, processes, the fillet rule, the flare default, the features-found check; the rc5 verification's deferred items (1)–(3); the E5 guard in `measure2`; the replica factories | the package |
| `apply_rc6_review_fixes.py` | the independent code review's three code items (§ 14.4): `devspecs` under `python -OO`, `features_found.by_kind`, the stretch-end refusal naming the outlet; comments | the package |
| `make_figure_anchors_r6.py` | `opencif/data/r6/figure_anchors.json` (hashed; the r5 file kept) | — |
| `apply_rc6_presets.py`, `apply_rc6_fig6_presets.py` | the presets and examples at the rc6 values | — |
| `measure_gifford_layout.py`, `apply_rc6_layout_presets.py` | rc6, after the author's review of 2026-10-17: the Gifford layouts the traces read on the figures without the mask (Fig. 4: the inset's leg spacing; Fig. 6: the leg axes, from the post rows and the turn walls' common centre, and the last leg's early end), gathered, then written into the five Gifford Fig. 4 and Fig. 6 presets | `gifford_layout_measurement.json` |
| `e5_default_check.py` | Reviewer A's E5 device (u05) at the new default fillet: every planned feature found; the build's own verdict (LIMITATIONS 23) | `e5_default_result.json` |
| `passage_scan.py` | the narrowest straight side channel of every shipped record against its cutoff D (LIMITATIONS 26; a report, opencif has no passage rule) | `passage_scan.json` |
| `dinh_length_bound.py`, `dinh_taper_alternatives.py`, `dinh_schedule_at_F.py` | rc5's scripts ported to take the preset's own pitch G + g_p and its tolerance (rc6: 45.1 ± 0.2 µm); `regress_length_rc5.py` shows that on the rc5 tree, at its 46.0 ± 1.0 µm, they reproduce rc5's outputs exactly | `dinh_length_bound_{depicted,specified}.json`, `dinh_taper_alternatives.json`, `dinh_schedule_at_F.json`, `out/*.out`, `out/length_regression_rc5.out` |
| `apply_dinh_notes.py` | the Dinh presets' `paper.note` and `paper.length_bound` written from those outputs (no number typed) | — |
| `compare_rc5_rc6.py`, `stats_rc6.py`, `band_rc6.py` | every record rc5 → rc6 (what moved and why), the record statistics the docs quote, the model band per design | `out/compare_rc5_rc6.*`, `out/stats_rc6.out`, `out/band_rc6.out` |

`rc6_traces/` holds the micrograph traces, re-run on rc6's masks, the whole-device figures and their code (its own
README); `RC6_DIFF.patch` is the diff from the rc5 post-review tree (10-12) to this one (records, logs, figures and PDFs
excluded).

## REVIEWS_2026-10-09.md

Summaries of the independent reviews of this round (code; code fixes; documentation; a third-party dry run on the
evaluation archive) with the disposition of every finding.  `docs/REVIEW_2026-10-08.md` holds the full history.
