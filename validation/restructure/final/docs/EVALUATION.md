# A guide for evaluating opencif 1.0.0

opencif generates Controlled Incremental Filtration (CIF) photomasks from two decoupled design targets, the size cutoff
D and the filtrate fraction F, and it checks every mask that it draws. The checks cover its own geometry, topology, and
design-rule checks on the drawn polygons, the one-dimensional flow network measured on them inside a stated validity
envelope, and the design code's physics auditor, ported, on the photomask file read back, with the auditor's own ruler
on the straight rows. These are statements of the model on the drawn mask, and they do not show that a chip was
fabricated or measured. This guide is written for an evaluator who decides whether the release can serve as the
dissertation's software, that is, as its Prospectus C7 (the targets-to-photomask tool) and C8 (the machine-enforced
validity envelope). It says what to run, what each result should be, where the evidence for every claim is, and what
is deliberately out of scope.

## 1. The materials

| location | contents |
|---|---|
| the repository | the package `opencif/`, the suites `tests/`, the design code `reference/CIF_v2.py`, `presets/`, `examples/`, `docs/`, `notebooks/`, `devices/`, `tools/`, `validation/`, and `provenance/` (`README.md` gives the map) |
| `validation/` | this release's summary of evidence, with `RELEASE.json` (the code and tests stamps and the counts), every photomask's sha256 (`photomasks_sha256.txt`) and RC8's (`rc8_comparison/frozen_photomasks.json`), the plans of RC8 and 1.0.0 and the record comparison (`rc8_comparison/`), every audit's result (`audit_table.csv`), the suites' logs, and the scripts that turned RC8's frozen tree into this one |
| the evidence assets of the release (`opencif-1.0.0_evidence_records.zip`, `opencif-1.0.0_evidence_figures_*.zip`, and `opencif-1.0.0_photomasks_*.zip`) | every build's record (`<name>.json`), summary, DRC result, and figures, the standalone audits, the build logs, and the photomasks, which are RC8's byte for byte (`validation/photomasks_sha256.txt`); they hold no publisher's figure. The device packages' outputs show the publishers' micrographs and are not distributed |
| `docs/OPTIONS.md` | the options sheet, which gives every `generate()` argument and YAML key with its code, default, origin, an example specification, and the record field that verifies it |
| `docs/PHYSICS_AUDIT.md`, `docs/CIF_V2_REQUIREMENTS.md` | the physics audit, check by check, and every requirement of the design code with its status in opencif |
| `docs/ALIGNMENT.md` | every statement of Appendix A, Appendix B, and the glossary on which the code relies, with its code, test, and status |
| `docs/LIMITATIONS.md` | the limitations and the future work |
| `docs/DISSERTATION_NOTES.md` | what the dissertation may claim, with the sources to cite |
| `docs/RESULTS.md`, `docs/results.csv` | every specification's verdict and key numbers, from one rebuild on one code stamp |
| `docs/REPRODUCTIONS.md`, `docs/PRESET_PROVENANCE.md` | each paper's printed quantities beside the reproduced ones, and the source of every preset value |
| `docs/REFERENCE.md`, `docs/ARCHITECTURE.md` | the full reference, and one build followed step by step |
| `docs/REVIEW_2026-10-08.md`, `docs/CHANGELOG.md` | the review history and the release history |

The texts that the code cites (Appendices A and B and the glossary, frozen on 2026-09-29, with the glossary's
preamble and the design specification v5) are in `foundation/`, and a copy of the preamble is in
`reference_glossary/`. `tests/test_appendix_trace.py` checks their md5 sums and confirms that the glossary values in
`reference/glossary_values.json` are the preamble's own.

## 2. Setup, and confirmation that the code is the code that made the evidence

```
python3 -m pip install -r requirements-tested.txt     # the exact versions with which every record was produced (Python 3.11.17)
python3 -m opencif --version                          # opencif 1.0.0 (code {{STAMP}})
python3 -c "from opencif import pipeline; print(pipeline.CODE, pipeline.environment_stamp()['tests_sha256'])"
```

The code stamp is a sha256 over every Python file of the package and its data, each hashed with its path relative to
the package (written with '/' on every system and taken in the same order), so any edit to the package changes it.
Every record of the evidence archive names it (`code.sha256`), and `validation/RELEASE.json` gives the tests' hash
({{TESTS_STAMP}}). The provenance scripts behind the package's two-dimensional tables also need
`provenance/rc5_scripts/requirements-provenance.txt`.

## 3. The tests

| suite | scope | checks | result | time |
|---|---|---|---|---|
{{SUITES_TABLE}}

The suites that need builds (`test_figure_anchors`, `test_rc7`, `test_rc8`, and `test_trace_anchors`) read `runs/`.
Rebuild first (section 4), or unpack the evidence archive's records and RC8's photomasks into `runs/`; without them,
these suites report INCOMPLETE (exit 4) rather than ALL OK. `test_trace_anchors` also needs the publishers' figures
(`OPENCIF_PAPER_RASTERS=notebooks/figures`, where the repository ships them; `provenance/rc6_traces/README.md`
explains how they were extracted),
and, for its taper-direction test (V13) of the three Fig. 4 variants, their builds: after a rebuild they are in `runs/`,
so set `OPENCIF_TRACE_VARIANTS=runs`.

Some checks have no seeded fault in the suites. They are the topology, side-channel continuity, and head-envelope
checks, the refusals for Re_gap alone and for the foot, an array refused for merging copies, a collision on the ports
alone, the stop conditions of the outlet extension, and a re-assessment that renames FAIL to PASS on a real record.
Those checks did fail on real drawings during development and caught real defects (REVIEW rows 38, 40 to 43, and 46),
and among the final records m14 is a real Q2 failure and `wbc_xia_process_refused` a real M13 refusal. Every hard and
envelope check of the physics auditor except P2 and P6b has a seeded fault in `tests/test_auditor.py`.

## 4. Rebuilding a device and comparing it with its record

```
python3 -m opencif build presets/gifford2014_fig4_design.yaml --out /tmp/check     # 3 to 10 min
python3 -m opencif build examples/matrix/m11_paperclip_const_pack_filter_ring.yaml --out /tmp/check
python3 -m opencif audit /tmp/check                                                # the audit again, from the file
```

Compare `/tmp/check/<name>.json` with the evidence archive's record. The verdict, `validation.realised_cutoff`,
`validation.network`, `validation.physics_audit`, and every `assessment.checks[*].value` should agree to the printed
precision with the tested versions. The DXF is written with fixed metadata, so with the tested versions the rebuilt
mask's sha256 equals the record's `files.dxf_sha256` (and RC8's). The command `python3 tools/run_batch.py 2
presets/*.yaml examples/*.yaml examples/*/*.yaml` rebuilds everything in about 5 h on two cores, and
`python3 tools/results.py` regenerates `docs/RESULTS.md`.

## 5. A new device from D and F

```python
import opencif
r = opencif.generate(3.0, 0.6, chip=(25, 75), layout='serpentine', inlet='filter', outlet='teardrop', out='/tmp/eval')
r.verdict, r.exit_code, r.physics_audit      # 'PASS ...', 0, 'PASS'   (a build takes 1 to 30 min)
r.dxf, r.figures, r.spec_path, r.record
```

The equivalent command is `python3 -m opencif generate 3 0.6 --chip 25x75 --layout serpentine --inlet filter
--outlet teardrop --out /tmp/eval`. D and F are the design targets, and every other argument is a part or the
fabrication context (`docs/OPTIONS.md`). With `check_only=True`, the call plans, gates, and draws the outline in memory
and writes nothing. The following calls show the refusals.

```python
fine = {'depth_um': 80, 'gap_um': 20.5, 'min_feature_um': 3, 'min_radius_um': 0.75}
opencif.generate(3.2, 0.8, width_in_um=100, meander=None, process=fine, check_only=True, out='/tmp/eval')
#   exit 5, refused on w_min (the straight side channel's 3.50 um at row 1, against h/AR_max = 8 um) and on the
#   kappa plateau (u_min 0.0320 < 0.035); band 2.777-3.200 um (App. B C7's method, gates.cutoff_band_um)
opencif.generate(4.0, 0.5, meander=None, process=fine, flow_uL_min=200, check_only=True, out='/tmp/eval')
#   exit 5, refused on w_min (4.10 um at row 1), on the channel Re (41.7 > 10) and on the Stokes gates
opencif.generate(4.0, 0.5, meander='hairpins', check_only=True, out='/tmp/eval')
#   exit 5, because the hairpin head cannot seal gap 1 here; the numbers are given, and no other head is substituted (about 1 min)
```

## 6. Claims and their evidence

| claim | evidence |
|---|---|
| D and F are the design targets (the channel's shape, the depth, and the gap and post sizes are the other inputs, whereas f_gap and N are outputs), and on the κ plateau they are independent (App. B thesis; App. A R13(c); App. B C8a) | `spec.SCHEMA` has no f_gap target; the traces "changing D changes N, not F", "changing F changes N, not the cutoff", and the N·D² invariance; `examples/decoupling/` in RESULTS.md; `tests/test_design_code.py` (the design code's M9, the schedule's independence of the flow, the viscosity, and the post length) |
| every gap of the drawn mask realizes D ± 1% and the device realizes F within 10⁻³ (Q2), and the drawing realizes its design (Q1), at nominal dimensions and nominal junction terms | every record's `assessment.checks` ("Q1 ...", "Q2 ..."); RESULTS.md; for the sensitivity, LIMITATIONS 13 to 15a and `validation.design_code_reports.V_fab_tolerance` |
| the photomask file passes the design code's physics auditor, ported, and the auditor's own ruler agrees with the generator's measurement | `validation.physics_audit` in every record; `python -m opencif audit`; `tests/test_auditor.py` (the port against CIF_v2's own auditor, and a seeded fault for every check); `docs/VALIDATION.md` (every shipped photomask, in the build and from scratch) |
| the physics is the frozen framework's (one kernel pair, the march, the closures, and the duct law), and the kernels are the design code's | `tests/test_appendix_trace.py`; `tests/test_kernels.py` (bit for bit); `docs/ALIGNMENT.md`; the records carry `design.kernel_pair` |
| every requirement of the design code is ported, reported, met by an equivalent, or not applicable, and the document says which | `docs/CIF_V2_REQUIREMENTS.md`; the trace "CIF_v2 Pairing / Tolerances / Process: every field accounted for"; `tests/test_design_code.py` |
| the validity envelope is enforced in strict mode (App. B thesis, with the plateau and the Reynolds and Stokes gates, together with M13, E11 on every row, M16, and the γ_GS and w_min gates checked before drawing); an outlet that cannot close the design is refused in every mode; and the chip area and the G9 leg-length budget are reported rather than enforced (App. B thesis; C8 CQ5) | `build.gates`; `pipeline.gamma_gs_check`, `pipeline.wmin_check`, and `pipeline.outlet_precheck`; the audit's M14 on the drawn mask; `test_negative.py`; `test_lattice.py`; the "envelope" column of RESULTS.md |
| below the plateau, the cutoff is reported as App. B C7's plateau-anchored band | `build.gates` (`cutoff_band_um`); the trace reproduces App. B's 1.25, 2.17, 4.08, and 1.80 µm |
| the papers' devices are reproduced from D and F (at the design level, rather than as a match to every printed number or measurement) | the presets in RESULTS.md; REPRODUCTIONS.md; PRESET_PROVENANCE.md; `test_figure_anchors.py`; the device kits (`devices/`) |
| the drawing has no intersections, the drawn parts are the ones requested, and no head is ever substituted | the records' checks; `test_negative.py`; `test_composition.py` (LIMITATIONS 29 lists what is not covered) |
| for fabrication, nothing narrower than the floor is lost, the fluid is one component, and in strict mode no CIF channel is deeper than `process.aspect_ratio_max` allows | "DRC at the process floor" in every record; the audit's P9 and P9-AR (LIMITATIONS 7) |
| the restructured package changes no design and no photomask | `validation/`, where every specification is planned by RC8 and by 1.0.0 (every array of the row schedule equal), every photomask is rebuilt and compared byte for byte, and every record is compared |
| the evidence comes from this code | one `code.sha256` across all records, and each mask's sha256 in its record |

## 7. What is deliberately not claimed

`docs/LIMITATIONS.md` lists every item. An evaluator should weigh the following ones first.

1. The realized cutoff is the one-dimensional network model on the drawn geometry rather than a measured particle
   cutoff, and no mask here has been fabricated.
2. The scope is rigid, approximately spherical particles in a Newtonian fluid of nearly constant viscosity (water in
   the Reynolds gates), with laminar channel flow that is creeping at the gap and lamina scales ("not globally
   Stokes", in the glossary's words). Platelets, deformable cells, crowding, and whole blood are out of scope (App. B
   thesis).
3. The design code's stricter inlet limit is reported rather than enforced. The physics auditor's P8 (the inlet Re
   against 1.0 and the naive entrance length against Λ/2) is checked on every photomask at the design flow, where every
   design exceeds the Re limit, and at CIF_v2's default of 5 µL/min, where every design passes. The design code's own
   gate M14 (Re ≤ 10) is enforced (LIMITATIONS 6), and the dissertation's erratum D106 says so.
4. The auditor's own ruler covers the straight rows, which hold most of every device's gaps. The head units, the
   U-turns, and the outlet rely on the generator's measurement, with their model terms and their lengths taken from the
   plan (LIMITATIONS 12).
5. One rule of spec v5, G11's fold of a side channel narrower than w_min, is replaced rather than implemented, and
   strict mode refuses such a design (LIMITATIONS 7).
6. The Dean (inertial) bend term is not modeled; a screening estimate applies where the bound on the Reynolds number of
   any channel is at most 10 (LIMITATIONS 11). The bend, head T-junction, and outlet stretch-end terms are
   two-dimensional fits, and the record's model band shows which rows their range moves.
7. The calibration interval of κ_g, ±3.4%, is wider than the design tolerance of ±1%.
8. D ± 1% is a statement about the nominal geometry. Per-gap cutoffs are very sensitive to width errors below a
   micrometer, as the records' scenarios (LIMITATIONS 13 and 14) and the design code's Monte Carlo fabrication tolerance
   (`V_fab_tolerance`, reported on every build) show.
9. No published measurement tests κ_g sharply, because Gifford's Fig. 5 devices retain about half or less of the beads
   of size D (`provenance/gifford_fig5_d50.py`).
10. The papers' own devices lie mostly outside the validity envelope, so they are reproduced in reproduce mode, with a
    cutoff band, and exit 3.
11. Gifford's Fig. 4 is reproduced on its claim (leg 4) and on its layout, whereas the neighboring leg's widths are not
    reproduced (LIMITATIONS 28).
12. The "inside the envelope" of `check` is not a prediction that the build passes Q2 (LIMITATIONS 30). Three known
    behaviors are left for a later release so that the photomasks of 1.0.0 stay RC8's. A head-less floor start outside
    reproduce mode is noted at the plan stage and fails Q2 at the build (m14), a plan whose F overshoots the tolerance
    is noted and fails Q2 unless `physics.trim_F` is set, and the leg fitter's collision search can be slow.

## 8. Suggested criteria for the freeze

1. Section 2 shows one code stamp across all records, equal to the package's, and the photomasks' sha256 sums equal
   RC8's.
2. Section 3 shows every suite ALL OK on a rebuilt tree (with the publishers' figures for `test_trace_anchors`).
3. Every photomask's physics audit is PASS or PASS_OUTSIDE_ENVELOPE, in the build and from scratch, apart from m14's
   designed FAIL (P5 on its undeclared floor-start rows), and PASS_OUTSIDE_ENVELOPE occurs exactly where the design is
   outside the validity envelope (`validation/audit_table.csv`).
4. RESULTS.md shows every preset REPRODUCED (exit 0) or REPRODUCED outside the envelope with a band (exit 3), and every
   new design exits 0 except the three deliberate outcomes (`dinh_targets_sealed_head`, exit 3;
   `m14_paperclip_const_nomeander_floorstart19`, FAIL, exit 2; and `wbc_xia_process_refused`, REFUSED, exit 5).
5. An independent re-audit of the rows of ALIGNMENT.md against the dissertation's texts finds no row whose status
   should be MISLABELLED or CONTRADICTED.
6. A spot rebuild (section 4) reproduces its record with the tested versions, and a new design from `generate()`
   (section 5) passes or is refused with reasons that match the gates.
