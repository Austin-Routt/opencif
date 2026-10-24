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
   scripts in `validation/restructure/` turn RC8's frozen tree into this one with every edit an exact match.
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
   `assessment.status.physics_auditor_P8` ("NOT RUN ON THE MASK") becomes `assessment.status.physics_audit`.
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
   and a seeded fault for every check), `test_teaching`, and `test_design_code` (each report against CIF_v2's own
   function, with the design code's determinism R1 and invariance M9 as tests). In `test_appendix_trace`, η and β_cal
   are now reported fields (each asserted equal to CIF_v2's), the disposition of L_eq names both of its replacements,
   and the glossary's numeric macros come from `reference/glossary_values.json` when the dissertation's texts are not
   present.
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
   revised for 1.0.0 (the dissertation's errata list; DISSERTATION_NOTES section 4).
10. The repository ships the texts that the code cites: Appendices A and B and the glossary (frozen on 2026-09-29, as
    LaTeX sources and as the PDFs compiled at the freeze), the glossary's preamble, the design specification v5, the
    errata list, and the differences from the previous freeze, in `foundation/` (`foundation/README.md`), so that
    `tests/test_appendix_trace.py` checks the texts themselves in every copy of the repository.
11. A plan that cannot be completed is a spec error (exit 6) with the design's numbers, rather than a traceback
    (exit 1). RC8 converted a closure without a solution in this way, whereas a march that needed more than
    400,000 rows (for example, D = 0.25 µm on a 100 µm center channel) stopped with the march's RuntimeError.
    `tests/test_negative.py` now seeds that case. The figures' titles separate the name from the targets with a
    colon. Neither change touches a shipped design or photomask.

Some things are unchanged by design, so that every photomask stays RC8's. They are the physics, the schedule, the
closure, the layout, and the drawing, together with three known behaviors left for a later release (LIMITATIONS 30).
A head-less floor start outside reproduce mode is noted at the plan stage and fails Q2 at the build, an F overshoot
beyond the tolerance is noted and fails Q2 unless `physics.trim_F` is set, and the leg fitter's collision search can be
slow. The British spellings of twelve specification and record keys (for example `centre_channel` and
`realised_cutoff`) and of the program's printed messages are also kept, because renaming the keys would break every
preset and record and the messages are stored in the records and compared by the tests; accepting American aliases on
input is left for a later release.

Validation: {{VALIDATION_SHORT}}
