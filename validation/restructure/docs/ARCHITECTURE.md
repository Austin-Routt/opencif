# One build, step by step

This page follows one build from the two targets to the audited photomask and names the module and the function that
carry out each step, so that a reader can open the code at the right place. `opencif/teaching.py` follows the same
path for the simplest device, a straight ladder, in fewer than 300 lines, and it is the place to start reading. The
full generator described below adds what real devices need, which is a meander head, tapers, folded layouts, inlets,
outlets, and the parts of the published devices.

```
targets D and F, process, parts
   │  opencif.spec               validate the specification and map it onto the engine's Spec
   ▼
plan                             opencif.physics.schedule, with opencif.physics.kernels
   │  δ* = D/κ_g;  f(i) = g(δ*/w_c(i));  march P and S until F;  the exact closure R_s(i);  the widths w_s(i)
   │  opencif.pipeline.prepare   the head's sealing search, the floor start, the validity envelope (strict mode refuses)
   ▼
draw the outline in memory       opencif.geometry.layout and opencif.geometry.shapes      (build.plan_and_draw)
   │  legs, head units, U-turns, inlet, and outlet; parts that collide are refused before anything is written
   ▼
write the photomask              opencif.geometry.dxf   (layer CHANNELS is the fluid and ISLANDS the posts; reproducible bytes)
   ▼
measure the drawn polygon        opencif.checks.measure  (cross-sections give the resistance of every drawn part)
   │  the drawn network gives each gap's emergent fraction f_i, its realized cutoff D_i = κ_g x(f_i) w_c(i), and F
   ▼
assess                           opencif.build, opencif.pipeline.assess, and opencif.checks.drc
   │  Q1, the mask draws its design; Q2, it realizes D and F; the DRC, topology, ports, clearances, corners,
   │  aspect ratio, and the validity envelope
   ▼
physics audit                    opencif.checks.auditor and opencif.checks.ruler  (the design code's auditor on the
   │                              file read back, with its own ruler on the straight rows)
   ▼
the design code's reports        opencif.checks.design_code  (recorded and never a pass condition)
   ▼
verdict, exit code, record (<name>.json), summary, and figures (opencif.figures.proof)
```

## 1. The specification (`opencif/spec.py`)

A YAML file names the targets, the process, an optional chip, and a choice of parts. The same specification can be
given as a dict, or written by `opencif.generate(D, F, ...)` in `opencif/api.py`. `spec.load` validates every key, so
that an unknown key or value is an error (exit 6), and `spec.to_engine` maps the specification onto the engine's
`Spec` (`opencif/physics/schedule.py`) and a dict of options (the gates' mode, the tolerances, and the parts' styles).
`docs/OPTIONS.md` lists every option with its default and its source.

## 2. The plan (`opencif/physics/schedule.py` and `opencif/physics/kernels.py`)

`kernels.py` holds the formulas of Appendix A on which the rest is built. They are the duct resistance r(w, h)
(Result 14), the skim cubic g(x) = 3x² − 2x³ and its inverse (Results 1 and 2), the cutoff pairing κ_g = 4.9138
(App. B, claim C3b), and the peak wall shear and the entrance length (spec v5). They are the design code's own
formulas, and `tests/test_kernels.py` shows that they return the numbers of `reference/CIF_v2.py` bit for bit.

`schedule.build` marches row by row. Each gap removes f(i) = g(δ*/w_c(i)) of the center flow on each side, so that
every gap has the cutoff D; the center flow after row i is P(i) = (1 − 2f) P(i−1); and the march stops at the first row
at which 1 − P reaches F (App. A, eq. march and eq. Frec). The function classifies the rows (the meander head's arms
and hairpins, the rectilinear transition, and the straight rows) and closes the design exactly, which means that each
side segment receives the resistance that carries its flow between two gaps, with the change of the gap drop counted
(spec v5, section 13). The exact closure reduces to App. A's PE closure when the gap drop does not change.
`schedule.plan` chooses the head's arms and the first gap's fraction. `pipeline.prepare` adds the policies around it
(the head's sealing search, `shrink_head`, and the floor start), the pre-flight checks (`gamma_gs_check`, `wmin_check`,
and `plateau_preflight`), and the validity envelope (`build.gates`, which evaluates the κ plateau, M13, M14, E11, and
M16). Strict mode refuses a design outside the envelope (exit 5) before anything is drawn.

## 3. The drawing (`opencif/geometry/`)

`pipeline.fit_layout` chooses the number of legs, which places a linear, serpentine, or paperclip layout on the chip.
`build.plan_and_draw` lays out the legs (`layout.plan_layout` and `layout.make_legs`), plans the head's units
(`shapes.plan_head`), solves the transition unit, and assembles the whole outline (`layout.assemble`). Every straight
side wall is drawn so that its integrated resistance equals the closure's. `build.collisions` checks the parts'
clearances on the outline in memory, so that parts that collide are refused (exit 5) and no photomask is written. The
published devices' depicted parts are in `geometry/depicted.py` and `geometry/dinh_inlet.py`.

## 4. The photomask (`opencif/geometry/dxf.py`)

`dxf.write_dxf` writes the fluid polygon, with its outline on layer CHANNELS, its holes (the posts) on layer ISLANDS,
lengths in micrometers, and every curve as a fine polyline. The file carries fixed metadata, so the same geometry
always gives the same bytes, and a photomask can therefore be pinned by its sha256 (`record.files.dxf_sha256`).

## 5. The measurement (`opencif/checks/measure.py`)

`measure.measure_device` cuts the drawn polygon at every part. It cuts the gaps at mid-wall, the side channels across
their whole length (integrating r(w)), the center channel at each post tip, the head units by area over width, the
U-turns along each side's own center line, and the outlet leads along their drawn paths. `measure.network` then solves
the three-rail network (the center channel and both sides, with no symmetry assumed) on those resistances and returns
each gap's emergent fraction f_i, its realized cutoff D_i = κ_g x(f_i) w_c(i), and the delivered F.

## 6. The assessment (`opencif/build.py` and `opencif/pipeline.py`)

`build.build` collects the validation blocks. They are the topology, the envelope, the corner radius, the clearances,
the continuity of the side channels, the drawn components against the request, the measured network, the intended
network (the closure together with the declared deliberate departures of a reproduction, from
`build.intended_and_declared`), the realized cutoff row by row, the model and fabrication scenarios
(`physics/scenarios.py`), the gates, the DRC of the written file (`checks/drc.py`, run in its own process), and the
physics audit. `pipeline.assess` turns the blocks into checks of four kinds. A hard check fails the build when it
fails (exit 2); an envelope check places the build outside the validity envelope (exit 3 in evaluate and reproduce
modes); an incomplete check gives exit 4; and the remaining checks are reported. Q1 asks whether the mask draws its
design, comparing the drawn network with the intended one at every gap, and Q2 asks whether it realizes the targets,
with every undeclared gap within D ± 1% and F within its tolerance.

## 7. The physics audit (`opencif/checks/auditor.py` and `opencif/checks/ruler.py`)

The auditor is the design code's own auditor (the `ladder_physics.audit` embedded in CIF_v2.py, with checks P1 to
P10), ported to opencif and run on the photomask file. The DXF is read back and must equal the drawn polygon vertex for
vertex (P0). The auditor's own ruler (`ruler.py`, which shares no code with `measure.py`) then measures every ordinary
straight row again from the file (the posts among the file's islands, the gap openings, the wall thickness, the row
pitch, and the center and side widths), and the generator's measurement must agree with it (R0). The auditor's
network uses the ruler's resistances on those rows and the generator's resistances for the head units, the U-turns,
and the outlet. The remaining checks cover the kernels' anchors (P1), mass conservation (P2), the closure realized by
the drawn geometry (P3), the network bulk against the kernel (P4), the uniformity of the lamina and the cutoff
guarantee (P5), the filtrate fraction and a closed form for it (P6 and P6b), port matching (P7), CIF_v2's own stricter
inlet limit (P8, reported), the design code's Reynolds and Stokes gates on the drawn mask (M14), the fabrication floors
(P9 and P9-AR), the posts (R6b, the design code's acceptance check, applied to every post), and the corner radii (P10).
`docs/PHYSICS_AUDIT.md` gives each check, its source in CIF_v2.py, the change in the port and its reason, and the
validation of the port. `python -m opencif audit RUN_DIR` repeats the audit of a finished build from scratch, because
it plans and draws the design again, reads the DXF, measures the file's polygon with both rulers, and runs the DRC on
the file, without loading any saved Python object.

The design code also attaches mandatory reports to every device and records operating validators in its acceptance
run. `opencif/checks/design_code.py` computes each of them by the design code's own formula on every build and writes
them to `validation.design_code_reports`. They are the pressure drop, ηF, the transfer bias, the sensitivity map, the
porosity, the operating regime at 5 µL/min and at the design flow, the steric screen, and a Monte Carlo fabrication
tolerance. As in the design code, none of them is a pass condition (`docs/CIF_V2_REQUIREMENTS.md`).

## 8. The verdict and the record

| exit code | verdict | meaning |
|---|---|---|
| 0 | PASS (or REPRODUCED) | every hard check passed, and the design is inside the validity envelope |
| 2 | FAIL | a hard check failed, and the photomask is kept as `<name>_FAILED.dxf` for diagnosis |
| 3 | PASS or REPRODUCED, outside the validity envelope | the hard checks passed, whereas a gate of the envelope did not (evaluate and reproduce modes) |
| 4 | INCOMPLETE | a check could not run (for example, the DRC was switched off) |
| 5 | REFUSED | the design was refused before drawing (strict mode outside the envelope, a head that cannot seal, an outlet that cannot close, or colliding parts) |
| 6 | SPEC ERROR | the specification is invalid |

The record `<name>.json` holds `input` (the specification as given), `design` (the plan), `validation` (every block
above, including `physics_audit`), `assessment`, `files` (the photomask with its sha256, and the figures), `code` (the
code stamp, a sha256 over every Python file of the package and its data), and `environment`. The `assessment` holds
`verdict`, `exit_code`, `checks`, `notes`, and `status`, which gives the verdict's parts as separate fields with
enumerated `status_codes`. `<name>_summary.md` is a readable digest of the record.

## 9. The module map

| module | lines | role |
|---|---|---|
| `opencif/__init__.py`, `api.py`, `__main__.py` | {{LINES:__init__.py,api.py,__main__.py}} | `generate(D, F, ...)` and `build_spec`; the command line (`build`, `check`, `audit`, and `generate`) |
| `spec.py` | {{LINES:spec.py}} | the specification, with its schema, validation, defaults, and mapping onto the engine |
| `pipeline.py` | {{LINES:pipeline.py}} | `prepare` (policies, pre-flight checks, and gates), `fit_layout`, `run`, `assess`, and the record with its summary |
| `build.py` | {{LINES:build.py}} | `plan_and_draw`, the validation blocks, the gates, the intended network, the declared rows, and `build` |
| `compose.py` | {{LINES:compose.py}} | parallel arrays of a validated element |
| `published.py` | {{LINES:published.py}} | the published devices as engine specifications |
| `teaching.py` | {{LINES:teaching.py}} | the whole path for a straight ladder, written to be read |
| `physics/kernels.py` | {{LINES:physics/kernels.py}} | Appendix A's kernels and the design code's supplied constants (equal to CIF_v2.py's, bit for bit) |
| `physics/schedule.py` | {{LINES:physics/schedule.py}} | the march, the head, the exact closure, and the bend and junction terms |
| `physics/network.py`, `scenarios.py` | {{LINES:physics/network.py,physics/scenarios.py}} | the symmetric network; the model and fabrication scenarios |
| `geometry/layout.py`, `shapes.py` | {{LINES:geometry/layout.py,geometry/shapes.py}} | the layout of the legs, head, U-turns, inlet, and outlet; the head's units |
| `geometry/depicted.py`, `dinh_inlet.py` | {{LINES:geometry/depicted.py,geometry/dinh_inlet.py}} | the published devices' parts as their figures draw them |
| `geometry/dxf.py` | {{LINES:geometry/dxf.py}} | the DXF writer |
| `checks/measure.py` | {{LINES:checks/measure.py}} | the measurement of the drawn polygon and its network |
| `checks/drc.py` | {{LINES:checks/drc.py}} | the design-rule check of a DXF (run as its own process) |
| `checks/auditor.py` | {{LINES:checks/auditor.py}} | the physics auditor (P0 to P10, R0, M14, and R6b) and the standalone audit |
| `checks/ruler.py` | {{LINES:checks/ruler.py}} | the auditor's own ruler, which measures the straight rows again from the file |
| `checks/design_code.py` | {{LINES:checks/design_code.py}} | the design code's mandatory reports and acceptance validators |
| `figures/proof.py`, `drawing.py` | {{LINES:figures/proof.py,figures/drawing.py}} | the three proof figures of each build |

Comments in the generator's modules keep the identifiers of the review items and checks that they answer (BUILD-1,
ADJ-C6, RP-B5, and the like), because the records and the documents cite them. The history of each change is in
`docs/CHANGELOG.md`.

## 10. The design code (`reference/CIF_v2.py`)

Appendix A names the design code `CIF_v2.py` (version 2.1.0). It is kept here byte for byte (md5
0e0961a50441afb94671b4470d41bace), and the package never imports it. The tests use it as the reference. The kernels are
proven equal to its own (`tests/test_kernels.py`), its physics auditor is run beside the port on its own demonstration
devices (`tests/test_auditor.py`), its reports are compared with the port's (`tests/test_design_code.py`), and the
results of Appendices A and B that it computes are pinned (`tests/test_appendix_trace.py`).
`docs/CIF_V2_REQUIREMENTS.md` lists every requirement that it checks and the requirement's status in opencif.
