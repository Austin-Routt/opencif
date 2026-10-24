# Validation of opencif 1.0.0

This page states what was checked for this release, how it was checked, and with what result, and `validation/`
holds the evidence. Everything here is a statement about the software and about the one-dimensional network model on
the drawn masks, because no device was fabricated or measured with this release (`docs/LIMITATIONS.md`).

The release has code stamp 9e4eb55a0affeef3 (a sha256 over every Python file of the package and its data, which every record
names) and tests stamp 6340f50956101e9d, and it was validated with Python 3.11.17 and the versions listed in
`requirements-tested.txt`.

## 1. The changes in 1.0.0 and the changes it must not make

opencif 1.0.0 is the eighth release candidate (RC8, code stamp 593d977fc83ec217) restructured into a package that an
undergraduate can read (`physics/`, `geometry/`, `checks/`, and `figures/`). The kernels are restated in
`physics/kernels.py`, the design code's physics auditor is ported and run on every photomask with its own ruler, a
teaching module and the `audit` command are added, a plan that cannot be completed is reported as a spec error, and the
comments and docstrings of every module are rewritten in the author's voice. The restructure is reproducible, because
the scripts in `validation/restructure/` turn RC8's frozen tree into this one through moves, import rewrites, the
auditor's integration, renamed words in comments, the rewritten comments, and the changes made after them
(`patch_final.py`), with every edit an exact match (`make_tree.sh` runs them in order). The documents' placeholders are
then filled from the evidence (`fill_docs.py`, with the values of `make_values_v1.py` and then, for the three passages
of sections 6 and 7 that wait for the device packages, the tutorial and the reviews, those of `make_values_late.py` in
`values_late.json`), the results table is written from the records (`tools/results.py`),
the notebooks are written by `make_tutorial_v1_voice.py` and `make_traces_v1.py`, and the device kits' manifests by
`make_device_manifests.py`. The tree that these steps make equals the published tree file for file, except that the
notebooks' cells carry new random IDs, the published tutorial is the executed one, and the six publishers' figures
in `notebooks/figures/` are added as they are (`NOTICE.md`). The rewrite of the comments is proven to change no code, because each module's abstract syntax tree, with its docstrings
removed, equals the tree before the rewrite, and every code token keeps its text and its column
(`validation/restructure/voice_check.py`). The restructure must change no design and no photomask, and the following
checks confirm that it changes none. RC8's photomask hashes, the plans of both releases, and the record comparison are
in `validation/rc8_comparison/`. The comparison covers the 53 records that RC8 also wrote (the three Fig. 4 trace
variants' photomasks and exit codes are compared, and their records are new). `validation/logs/rebuild_summary.jsonl`
also holds the 25 rows of the first rebuild pass, at code 8e9c4de5c2300b3f, before the ruler's corners were corrected
(CHANGELOG item 12); the last row of each design is the result.

| check | result |
|---|---|
| every shipped specification planned by RC8 and by 1.0.0, comparing the specification, every array of the row schedule (by sha256), the gates, the options, and the notes (54 specifications) | 54 of 54 identical in every compared field; the two fields compared apart are the planning time and the text of gate P8's status, which 1.0.0 rewrites on all 54 because it now audits every photomask |
| every shipped specification built again with 1.0.0 (54 builds, 53 photomasks, and 2 arrays), with each photomask compared with RC8's byte for byte | 53 of 53 photomasks and 2 of 2 arrays identical to RC8's, byte for byte (sha256 in `validation/photomasks_sha256.txt`); the refused design (`wbc_xia_process_refused`) draws no photomask in either release |
| every record compared with RC8's, except for what 1.0.0 adds or renames, each named with its CHANGELOG item (the physics audit and its status fields, the design code's reports, an array's element audit, the Pairing's new module, the text of gate P8's status, and the module names in two notes), and the timings | 52 of 53 identical apart from the named additions; 1 changed for a named reason (`m14_paperclip_const_nomeander_floorstart19`: its validity envelope goes from inside to outside, because the audit's M14 on the drawn mask finds the floor start's undeclared rows above the Stokes gate, CHANGELOG item 3; its exit code, 2, is unchanged); 0 unexplained (`validation/RELEASE.json`, counts.records) |
| the exit code of every build | 54 of 54 as in RC8 (exit 0: 40, exit 2: 1, exit 3: 12, exit 5: 1) |

## 2. The physics audit of every photomask

| audit | result |
|---|---|
| in the build (53 photomasks) | 47 PASS, 5 PASS outside the envelope (`gifford2014_fig5D_as_fabricated_12p5um`, `xia2016_rs7_aligned`, `xia2016_rs7_depicted`, `xia2016_rs7_depicted_pullR`, `xia2016_rs7_specified`: M14 on the drawn mask), and 1 FAIL (`m14_paperclip_const_nomeander_floorstart19`: hard P5, envelope M14, the design that already fails Q2); the refused design has no photomask to audit |
| again from scratch (`opencif audit`, which plans and draws the design again, reads the DXF back, measures the file's polygon with both rulers, and runs the DRC on the file) | the same status as in the build for 53 of 53 photomasks (47 PASS, 5 PASS outside the envelope, 1 FAIL) |
| the agreement of the two rulers on the straight rows (R0) | a relative difference of at most 1.3 × 10⁻⁵ in the segment resistances and 2.0 × 10⁻⁸ µm in the gap positions, over 132,991 straight segments and 266,328 gaps |
| the largest values over the strict designs | over the 27 strict designs: the drawn fraction of an undeclared gap within 0.47% of the intended one (P4, `dc_D3_F0p5`); the lamina δ/δ* from 0.9982 to 1.0023 (P5); the exact closure within 4.6 × 10⁻⁵ (P3); the drawn F within 4.0 × 10⁻⁴ of its target (P6); every post within 1.2 × 10⁻¹¹ µm of its planned size (R6b); and the lamina Reynolds number at most 8.8 × 10⁻³ (M14, against 10⁻²) |
| the reproductions | 18 presets and the 3 Fig. 4 variants: 16 PASS, 5 PASS outside the envelope (`gifford2014_fig5D_as_fabricated_12p5um`, `xia2016_rs7_aligned`, `xia2016_rs7_depicted`, `xia2016_rs7_depicted_pullR`, `xia2016_rs7_specified`); their declared departures are reported row by row, and a reproduction outside the envelope is the paper's design, not opencif's choice (`docs/REPRODUCTIONS.md`) |

`validation/audit_table.csv` gives every check of every photomask.

## 3. The test suites

The suites ran at code 9e4eb55a0affeef3 and tests stamp 6340f50956101e9d, except that the builds' records carry the tests stamp
5d83c0644e5bb4e3, written before `tests/test_composition.py` gained its audit allowance (CHANGELOG item 6). Every log
is in `validation/logs/`, with the tree's own folder written `<repo>`.

| suite | scope | checks | result | time |
|---|---|---|---|---|
| `test_kernels` | every kernel equals the design code's, bit for bit, over grids of inputs | 18 | ALL OK | 6 s |
| `test_teaching` | the teaching module against the kernels | 13 | ALL OK | 4 s |
| `test_units` | every unit conversion and dimensionless group of the record, and the curved-duct limits | 34 | ALL OK | 154 s |
| `test_bend_fd` | the curved-duct closed form against a finite-difference solve | 3 | ALL OK | 6 s |
| `test_design_code` | the design code's reports and validators, each against CIF_v2's own function | 18 | ALL OK | 126 s |
| `test_auditor` | the ported auditor against CIF_v2's, its solver, fourteen seeded faults, the ruler's corners | 51 | ALL OK | 369 s |
| `test_negative` | refusals and spec errors: every bad input caught with its exit code | 64 | ALL OK | 3051 s |
| `test_lattice` | the plan over a lattice of targets and processes | 42 | ALL OK | 484 s |
| `test_presets_match_replicas` | each preset plans its replica's design | 5 | ALL OK | 943 s |
| `test_composition` (two halves) | the parts composed over 22 cases and extras | 53 + 43 | ALL OK | 1614 s + 3086 s |
| `test_rc7` | the rc7 acceptance checks on the photomasks and trace results | 112 | ALL OK | 271 s |
| `test_rc8` | the RC8 acceptance checks, with the arrays recomposed | 67 | ALL OK | 2191 s |
| `test_figure_anchors` | the reproductions against the papers' printed figures (needs the papers' page images) | 42 | ALL OK | 94 s |
| `test_trace_anchors` | the traces against the micrographs (needs the papers' page images and the variant builds) | 339 | ALL OK | 2632 s |
| `test_options_sheet` | docs/OPTIONS.md against the specification's keys | 20 | ALL OK | 2 s |
| `test_appendix_trace` | every constant and identity the code uses against the shipped derivation documents | 73 | ALL OK | 13 s |

## 4. The auditor

`tests/test_auditor.py` tests the port against CIF_v2.py's own auditor on its three demonstration devices, which give
the same numbers and the same verdicts, and tests its solver against the build's solver on random ladders. It then
seeds fourteen faults, one at a time, into a drawn device and confirms that each is caught by the check meant to catch
it (`docs/PHYSICS_AUDIT.md`). Last, it checks the ruler's corners on a synthetic taper that changes slope at one gap.

## 5. The kernels

`tests/test_kernels.py` shows that every kernel equals the design code's, bit for bit, over grids of inputs. The grids
cover the duct bracket at 3,001 aspect ratios, the duct resistance at 2,400 pairs of width and depth, its inverse at
300 resistances, the cubic at 20,001 points, its inverse at 4,004 fractions, the shear coefficient, the entrance
length, and the post outline. The suite also shows that every field of the pairing, tolerance, and process records is
equal and that the design code is the frozen v2.1.0 file.

## 6. The tutorial and the device packages

The tutorial (`notebooks/opencif_tutorial.ipynb`) was executed in full (`QUICK = False`) on this tree with the tested versions in 50 min, and none of its 37 code cells is in error. The published notebook is that executed one, with the text of its section 12 updated afterwards to say that the figures are distributed (no code cell or output changed), and it contains no absolute path of the machine that ran it. Its section 6 built a new device from D = 4 µm and F = 0.5 (verdict PASS, exit 0, physics audit PASS). Its section 7.1 audited that photomask again from scratch, and the audit gave PASS, with every check's verdict equal to that of the build's own audit (the checks whose recorded value differs in its last digits or in the row where it peaks: R0, P2, P4, P5). Its section 10 rebuilt the four published devices, and every photomask is byte-identical to the release's.

The four device packages were assembled from this tree by `tools/device_package.py` and run with the publishers' figures on 2026-10-23, each in its own process. Each first checked the tested versions, the code stamp and the sha256 of every input, then built its photomasks, registered them on the papers' micrographs and rendered the overlays, and every one ended its last step. All 11 photomasks are byte-identical to the release's, which are RC8's. Per package: `gifford2014_fig6` (43 min): two photomasks byte-identical, four registration records equal to RC8's in every value, 20 measurements within their regression limits; `gifford2014_fig4` (37 min): five photomasks byte-identical, five registration records equal to RC8's in every value, every result of the taper-direction test (V13) equal to RC8's, with the nested reading's neighbouring leg PASS, reading B's FAIL and every leg-4 control PASS, the claim's leg on the inset within its regression limits; `dinh2024_v1` (21 min): two photomasks byte-identical, two registration records equal to RC8's in every value; `xia2016_rs7` (54 min): two photomasks byte-identical, 12 registration records equal to RC8's in every value. The packages' outputs, which embed the publishers' figures, are among the release's evidence assets rather than in this repository. The first run of the Xia package stopped at its trace step on an error of the restructure in five trace scripts, which was corrected (CHANGELOG item 14) before all four packages were run again; the results above are from that second run. After it, the kits' source notes (`inputs/SOURCE.md`) were updated to say that the figures are distributed in `notebooks/figures/`, which changed the sha256 of each kit's source note in its `MANIFEST.json` and `reference/EXPECT.json` and no other entry, and each kit then assembled from the published tree with the figures of `notebooks/figures/` and passed its check step.

## 7. Independent review

Two reviewers worked independently on the tree as published and on the build's evidence (2026-10-23): one on the physics, the drawing and its measurement, and the gates, and the other on reproducibility, the documents against the records, the tutorial, and the repository's hygiene. Both found the physics sound. Every finding was checked against the code or the records before it was acted on. The release keeps its code stamp, because its photomasks, records and audits were built and validated at that stamp, so a finding about the package is answered by documents that say what the code does, and the code change is listed for the next release (`docs/LIMITATIONS.md` item 31). The documents, the tests' workflow, the tools outside the package and the provenance were corrected where needed (CHANGELOG item 14). A third reviewer re-audited every row of `docs/ALIGNMENT.md` against `foundation/` and the restructured code (criterion 5 of `docs/EVALUATION.md` section 8). It found one row to relabel, AA16, whose floor start departed from spec v5 §10 and M6 (row S3, erratum D108), and ten smaller corrections; after they were made, its re-check found every row correctly labelled. At the author's direction the specification was then corrected in place (revision 5.1, `foundation/README.md`). The two suites that read the documents were run again after every change, with identical logs. The reviews and the record of what was done with each finding are among the release's evidence assets. With these, the release meets the six criteria of `docs/EVALUATION.md` section 8: sections 1 to 3 of this page and `docs/RESULTS.md` show criteria 1 to 4, the re-audit shows criterion 5, and for criterion 6 the commands of its sections 4 and 5 were run as printed: the rebuilt `gifford2014_fig4_design` and `m11_paperclip_const_pack_filter_ring` equal their records in every compared field (the verdict, the realised cutoff, the network, the physics audit, every check's value, and the photomask's sha256), and the new design passes and the three refusals refuse for the reasons printed (the comments on two of them were corrected to name the refusal on w_min at row 1 that comes first).
