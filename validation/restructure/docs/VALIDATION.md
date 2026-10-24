# Validation of opencif 1.0.0

This page states what was checked for this release, how it was checked, and with what result, and `validation/`
holds the evidence. Everything here is a statement about the software and about the one-dimensional network model on
the drawn masks, because no device was fabricated or measured with this release (`docs/LIMITATIONS.md`).

The release has code stamp {{STAMP}} (a sha256 over every Python file of the package and its data, which every record
names) and tests stamp {{TESTS_STAMP}}, and it was validated with Python {{PYTHON}} and the versions listed in
`requirements-tested.txt`.

## 1. The changes in 1.0.0 and the changes it must not make

opencif 1.0.0 is the eighth release candidate (RC8, code stamp 593d977fc83ec217) restructured into a package that an
undergraduate can read (`physics/`, `geometry/`, `checks/`, and `figures/`). The kernels are restated in
`physics/kernels.py`, the design code's physics auditor is ported and run on every photomask with its own ruler, a
teaching module and the `audit` command are added, a plan that cannot be completed is reported as a spec error, and the
comments and docstrings of every module are rewritten in the author's voice. The restructure is reproducible, because
the scripts in `validation/restructure/` turn RC8's frozen tree into this one through moves, import rewrites, the
auditor's integration, renamed words in comments, and the rewritten comments, with every edit an exact match. The
rewrite of the comments is proven to change no code, because each module's abstract syntax tree, with its docstrings
removed, equals the tree before the rewrite, and every code token keeps its text and its column
(`validation/restructure/voice_check.py`). The restructure must change no design and no photomask, and the following
checks confirm that it changes none.

| check | result |
|---|---|
| every shipped specification planned by RC8 and by 1.0.0, comparing the specification, every array of the row schedule (by sha256), the gates, the options, and the notes ({{N_SPECS}} specifications) | {{PLAN_RESULT}} |
| every shipped specification built again with 1.0.0 ({{N_BUILDS}} builds, {{N_MASKS}} photomasks, and {{N_ARRAYS}} arrays), with each photomask compared with RC8's byte for byte | {{MASK_RESULT}} |
| every record compared with RC8's, except for the new physics audit, the renamed status field, the module names in two notes, and the timings | {{RECORD_RESULT}} |
| the exit code of every build | {{EXIT_RESULT}} |

## 2. The physics audit of every photomask

| audit | result |
|---|---|
| in the build ({{N_MASKS}} photomasks) | {{AUDIT_BUILD}} |
| again from scratch (`opencif audit`, which plans and draws the design again, reads the DXF back, measures the file's polygon with both rulers, and runs the DRC on the file) | {{AUDIT_STANDALONE}} |
| the agreement of the two rulers on the straight rows (R0) | {{R0_SUMMARY}} |
| the largest values over the strict designs | {{AUDIT_STRICT}} |
| the reproductions | {{AUDIT_REPRO}} |

`validation/audit_table.csv` gives every check of every photomask.

## 3. The test suites

| suite | scope | checks | result | time |
|---|---|---|---|---|
{{SUITES_TABLE}}

## 4. The auditor

`tests/test_auditor.py` tests the port against CIF_v2.py's own auditor on its three demonstration devices, which give
the same numbers and the same verdicts, and tests its solver against the build's solver on random ladders. It then
seeds fourteen faults, one at a time, into a drawn device and confirms that each is caught by the check meant to catch
it (`docs/PHYSICS_AUDIT.md`).

## 5. The kernels

`tests/test_kernels.py` shows that every kernel equals the design code's, bit for bit, over grids of inputs. The grids
cover the duct bracket at 3,001 aspect ratios, the duct resistance at 2,400 pairs of width and depth, its inverse at
300 resistances, the cubic at 20,001 points, its inverse at 4,004 fractions, the shear coefficient, the entrance
length, and the post outline. The suite also shows that every field of the pairing, tolerance, and process records is
equal and that the design code is the frozen v2.1.0 file.

## 6. The tutorial and the device packages

{{TUTORIAL_RESULT}}

{{DEVICES_RESULT}}

## 7. Independent review

{{REVIEW_RESULT}}
