# The physics audit of every photomask

Every photomask that opencif writes is read back from its file and audited. The result is
`validation.physics_audit` in each build's record, and `python -m opencif audit RUN_DIR` repeats the audit for a
finished build. The auditor (`opencif/checks/auditor.py`) is the design code's own physics auditor, ported to read a
photomask file, with its own ruler for the straight rows (`opencif/checks/ruler.py`).

## Why the audit exists

Appendix B relies on the design code's physics auditor twice. Its scope discussion states that "the physics auditor
separately checks the drawn device's inlet against its own, stricter channel limit" (App. B, lines 4349 and 4350), and
it states that the targets-to-photomask tool and the machine-enforced validity envelope (Prospectus C7 and C8) "are
delivered as verified software, warranted by their acceptance suite and the physics auditor rather than by a
defeasible inference" (lines 4375 to 4377). The auditor is the `audit()` function of `ladder_physics`, a module that
the design code `reference/CIF_v2.py` (v2.1.0) embeds as `_PHYSICS_SRC`, with ten checks, P1 to P10.

No release before 1.0.0 ran it. CIF_v2's auditor reads CIF_v2's own meander patterns (the stroke widths that they
declare, their trunk outline, and their SVG fillets), so it cannot read a DXF file or opencif's parts, and its P8 limit
(an inlet Reynolds number of at most 1 at the device's flow) would refuse every opencif device at its real flow. From
rc4 onward the records said so (rc4 recorded "ABSENT", and rc5 to RC8 recorded "NOT RUN ON THE MASK" in
`assessment.status.physics_auditor_P8`), and open item D95, recorded in the rc4 review round, proposed naming CIF_v2 in
App. B's sentence. opencif 1.0.0 ports the auditor and runs it on the photomask file of every build. Because the port
enforces the design code's gate M14 on the drawn mask and reports the auditor's stricter P8 limit rather than enforcing
it, the dissertation's erratum names both. Thus, D95 is closed into D106, revised for 1.0.0, and
`docs/DISSERTATION_NOTES.md` (section 4) gives the wording.

## What the auditor reads

The auditor reads four things, in this order.

1. The file. The DXF is read back (layer CHANNELS holds the outlines of the fluid, and layer ISLANDS holds the posts
   inside them) and must equal the drawn polygon vertex for vertex, in micrometers, with closed polylines of straight
   segments only (P0). Everything measured afterward is therefore a measurement of the file.
2. Two rulers on the straight rows. The generator measures its own drawing (`opencif.checks.measure`), reading every
   part's widths where the plan puts it and taking the parts' lengths from the plan. An auditor that used only that
   measurement would share any error in it, so the auditor measures every ordinary straight row again with its own
   ruler (`opencif.checks.ruler`), which shares no code with the generator's measurement. The ruler finds every post
   among the file's islands and reads each gap's opening on five cross-sections across the wall, the wall's thickness
   from the posts' extent across the flow, the row pitch from one gap's center to the next, the center width at every
   post tip, and the side width from the posts' outer tips to the outer wall, integrated piecewise between the
   polygon's own vertices. The plan serves only to identify the ordinary straight rows and to say where to start
   looking for each post. The two measurements must agree (R0), and the auditor's network uses the ruler's
   resistances on those rows.
3. The other parts. The head units, the U-turns, and the outlet are taken from the generator's measurement. They carry
   the model's own terms (the curved-duct correction, the junction terms, and the outlet's end term), and their lengths
   and arc radii come from the plan (`docs/LIMITATIONS.md`, item 12), so the auditor's ruler does not measure them.
4. Its own network. The auditor assembles the resistor network of the drawn device itself (`auditor.solve_ladder`,
   with its own node numbering, solved by a banded Cholesky factorization) rather than using the build's solver
   (`measure.solve3`, a sparse LU factorization), and it compares the two on the same resistances (P2).

The standalone audit (`opencif audit RUN_DIR`, which calls `auditor.audit_run`) starts from the record's input alone.
It plans and draws the design again in memory, reads the DXF, requires the file to equal that drawing (P0), measures
the polygon read from the file with both rulers, runs the DRC on the file, and loads no saved Python object.

## The checks

Each check has one of three roles. A failure of a hard check fails the build (exit 2) in every mode. A failure of an
envelope check puts the drawn device outside the validity envelope (exit 3), and strict mode refuses such designs
before drawing them. A report is computed and recorded and is never a pass condition. The declared departures (a
reproduction's depicted head, a paper's transition, a floor start, and an allowed seed) are the rows that
`build.declared_rows` names from the design before the drawing is measured. P3 to P5 and the gap terms of M14 are
evaluated on the other rows, and the declared rows are reported apart.

| check | requirement | role | origin in CIF_v2 | change in the port, and the reason |
|---|---|---|---|---|
| P0 | the DXF read back equals the drawn polygon, vertex for vertex | hard | none | new, because CIF_v2 audits a pattern in memory and has no file to read |
| R0 | every post of every ordinary straight row is found on the file, and every straight row's center segment, side segment, and gap, as measured by the auditor's own ruler, is within 10⁻³ of the generator's measurement | hard | none | new, so that the auditor does not depend on the code that it audits (see above) |
| P1 | kernel anchors B(1) = 0.4217310, a square duct's R w h³/(μL) = 12/B(1) = 28.454154, g(0.3) = 0.216, g(½) = ½, g(x(f)) = f, and κ_g = 4.9138 | hard | `P1_kernel_anchors`, with B(1), c₂(1), the cubic limit of its laminar kernel, and the hard bound | the anchors of the kernel pair that opencif uses, the cubic with κ_g (App. A's pairing rule, Result 9; App. B, claim C3b), whereas CIF_v2's physics module uses its laminar kernel F_lam with κ = 5.798 |
| P2 | the plan's mass identity \|P(N) + 2[S(N) − σ₀] − 1\| < 10⁻¹⁰ (App. A); on the drawn network, the current law at every node to 10⁻⁴ of the smallest gap flow at an undeclared row, and the three outlets carrying the inlet flow to that tolerance times the 3N nodes, because their sum differs from the inlet flow by exactly the sum of the node residuals; 1 − F = ∏(1 − f_A − f_B) (App. A, eq. Frec) to 10⁻⁹; and, on the generator's resistances, the build's network solution within 10⁻⁶ (of the larger of f and f*) of the auditor's | hard | `P2_mass_identity`, the schedule's identity | the drawn network's own conservation is added, with a tolerance tied to the smallest flow that it measures (a network of 24,000 nodes keeps about 10⁻⁹ in double precision, which is why the outlets' tolerance grows with the number of nodes), together with the agreement of two independent solvers |
| P3 | every drawn side resistance is within 0.15 of the design's exact closure, recomputed from the drawn center and gap resistances, and of the PE closure (P/S) R_c (App. A, eq. closure) | hard | `P3_closure_drawn`, which compares R_exc with ρ R_c within a sanity bound of 0.15 | the same bound, applied to opencif's own closure (the exact gap-drop closure of spec v5, section 13) as well as to the PE form that CIF_v2 checks, on the segments between undeclared gaps and on the outlets |
| P4 | \|f_i / f_i,int − 1\| ≤ 5 × 10⁻³ over the bulk, where f_int is the intended network's fraction, which is the kernel g(δ*/w_c(i)) at every row of a design without declared departures; the bulk is every row except the first, the terminal window n_term = max(2, ⌈λ ln(1/ε_net)⌉ + 1), and the rows beside each declared departure that its disturbance needs to decay, max(n_term, ⌈λ ln(A/ε_net)⌉ + 1) on each side; here λ = 1/ln z₊, z₊ is the root above 1 of z + 1/z = 2 + q, q = min (2R_c + R_s)/R_g, ε_net = 10⁻³, and A is the departure's largest deviation; the deviation from the kernel itself is reported (`kernel_bulk_max`) | hard | `P4_network_bulk`, with the same tolerance, window, and bulk, against the kernel | two changes for the declared departures, which CIF_v2's devices do not have. First, a window on both sides of every declared departure, because a departure disturbs the network as a boundary does, and a departure whose fractions are off by A needs ln A/ln z₊ more rows than CIF_v2's window. Second, the reference is the intended network (`build.intended_network`, the closure with the declared departures as intended), because a departure can shift the rows downstream of it without decaying. On Xia's reproductions the loop pull-in narrows the side channels of every row after the first U-turn, so that the undeclared rows of the third leg carry fractions about 2% below the kernel (cutoffs about 1% below D, within Q2) and those of the later legs up to 1% above it. The intended network carries the same shift, and against it the drawn network of Xia's aligned reproduction agrees to 6.4 × 10⁻⁶ over the bulk (and to 3.1 × 10⁻⁴ on every undeclared row, the windows included), whereas the kernel deviation is 1.9%. For a design without declared departures the reference is the kernel, as in CIF_v2 |
| P5 | δ_i/δ* = x(f_i) w_c(i)/δ* ≤ 1.02 at every undeclared gap on both sides, each with 0 < f < ½ (the cutoff guarantee) | hard | `P5_lamina_uniform`, with the same tolerance | both sides of the two-sided network, with the declared rows reported apart |
| P6 | F_drawn ≥ F_target − 10⁻³ and \|F_drawn − F_march\| < 0.02 | hard | `P6_F_target`, the same | none |
| P6b | a closed form for F, App. A's march-outlet identity F = 2R_c/(R_s + 2R_c) at the drawn outlet (per side), within max(8 × 10⁻³, 4 f̄) of the drawn network's F, with the identity's premise (PE at the outlet row) and its residual reported | hard | `P6b_F_closed_corrected`, its corrected continuum closed form for its own canonical taper (a ratio of c₂-weighted widths), with that tolerance | CIF_v2's closed form exists only for its canonical laminar taper, whereas the identity that App. A proves "for every gap schedule" (the result that the march and the outlet identity give the same F) applies to every opencif design; the tolerance is CIF_v2's |
| P7 | the F and the largest cutoff that result when both side outlets' resistance is tripled | report | `P7_port_matching` (the rail run-out tripled), informational | both sides of the two-sided network |
| P8 | CIF_v2's P8, the auditor's own and stricter inlet limit, which requires the drawn inlet's Re ≤ 1 and the naive entrance length 0.05 Re D_h ≤ Λ/2, evaluated at the design flow and at 5 µL/min | report | `P8_stokes_gates`, the same two criteria at its `Q_in_uLmin` (5 µL/min by default), a pass condition there | App. B states that the auditor "checks ... against its own, stricter channel limit", so the port checks that limit on the drawn inlet under its own name and records the result. It is not a pass condition, because App. B gives the design code the channel gate Re ≤ 10 (M14, the next row); the operating points that App. B discusses (Dinh's V1 at 30 µL/min, with a channel Re of about 5) lie inside that gate and outside Re ≤ 1; and 0.05 Re D_h is the entrance length that spec v5, the glossary, and CIF_v2's own `entrance_length` call wrong in this regime (`docs/LIMITATIONS.md`, item 6) |
| M14 | the design code's Reynolds and Stokes gates on the drawn mask at the design flow, which require the center channel's Re ≤ 10 at every row and Re_gap and Re_δ ≤ 10⁻² at every cutoff gap, computed from the drawn widths and the drawn network's flows | envelope | the gate M14 of `design()`, applied to the design | the gate that App. B gives the design code, evaluated on the drawn mask as well as on the plan |
| P9 | the DRC at the process floor finds no lost feature and one fluid component, and the narrowest drawn CIF channel is at least F_min wide | hard | `P9_fab_floors`, which requires the narrowest stroke ≥ w_min, the wall ≥ wall_min, and the minimum clearance | measured on the polygon, because a DXF has no strokes, so the DRC tests every open and solid feature at the floor |
| P9-AR | the depth-to-width ratio of the narrowest drawn CIF channel is at most AR_max (CIF_v2's w_min = max(F_min, h/AR_max)) | hard in strict mode, a report otherwise | part of `P9` (w_min) | as in `pipeline.assess`, a pass condition in strict mode and a report in evaluate and reproduce modes (`docs/LIMITATIONS.md`, item 7) |
| R6b | every post of every straight row is read back from the file (the head's rows and the rows next to a U-turn are left out), and its extent across the flow is within 0.5 µm of the wall thickness ℓ_g, its extent along the flow is within 0.5 µm of g_p, and its ratio of the two is within 0.05 of ℓ_g/g_p | hard | the acceptance validator `R6b` ("the DXF posts match the original Gifford 2:1 pill"), with the same three criteria on one post of its own DXF | every post rather than one; it also measures on the file the wall thickness ℓ_g that spec v5's gap resistance takes as the gap's length |
| P10 | the smallest local radius of curvature of the drawn boundary is at least the specification's r_min | hard | `P10_rounded_corners`, which requires round joins and caps to be declared, one fillet per concave corner to be rendered, and the corner radius to be within its cap | measured on the polygon's boundary, because a DXF has no fillet primitives to count |

The tolerances that P1 to P10 share with CIF_v2 are CIF_v2's own. `tests/test_auditor.py` reads them from the frozen
source text and requires the port's values to be equal.

## Validation of the port

The port was validated in five ways.

First, against the original (`tests/test_auditor.py`, part 1), CIF_v2's own auditor is run on its three demonstration
devices, with 69, 88, and 3,139 units. Given the same measured geometry and the same kernel, the port's functions
return the same terminal window and bulk, the same P3 (its PE closure), P4, P5 (the lamina ratio at every unit), P7,
and P8 values, and the same verdicts, including P6's. The port's network solver reproduces CIF_v2's network to the
precision that either solver reaches on it. Against an 80-bit extended-precision solution of the same network, the
port is within 3.7 × 10⁻⁸ and CIF_v2's own solver within 4.0 × 10⁻⁸ on the largest device, because CIF_v2's network
floats at the pressure that its 1-ohm ground tie gives the outlet, so double precision keeps about 10⁻⁸ of a gap's
flow. P1, P2, P6b, P9, R6b, and P10 use other formulas (the table above) and are validated by parts 2 and 3.

Second, against the build's solver (part 2), the auditor's banded Cholesky solution and `measure.solve3` agree to
4.5 × 10⁻¹⁰ on random ladders, and to {{PAR_DEVICES}} on every shipped device (P2).

Third, the two rulers (R0) agree on the shipped devices to {{R0_SUMMARY}}.

Fourth, seeded faults (part 3) show that each check catches the fault that it is meant to catch. A drawn device first
passes every check, and then one fault at a time is introduced. An intended network 1% away from the drawn one fails P4
alone, which shows that P4 compares the drawn fractions with the reference that it is given. Moving one vertex of the
file by 1 nm fails P0. Shifting the cubic kernel by 10⁻⁹ fails P1. Widening one gap by 0.6 µm fails P4 and P5 at that
row, R6b (the post beside it becomes 0.6 µm shorter), and R0, because the gap's center moves by 0.3 µm, which the
ruler's pitch detects and the generator's lengths, taken from the plan, do not. Narrowing one side segment to 60% fails
P3 at that segment, together with P4 and P5, and M14 places the device outside the envelope. A post 0.6 µm thicker
across the flow fails R6b at that post, together with P4, P5, and R0. A post moved half a gap along the flow fails R0,
because the generator takes the row pitch from the plan whereas the ruler takes it from the posts, and also fails P4 and
P5. A side segment made 2% narrower over half its length fails P4, because its gap's fraction moves by 47%, which is the
sensitivity of the rows that `docs/LIMITATIONS.md` (item 13) describes. It also fails P5, and it fails R0 because the
generator samples the stepped wall at fixed points whereas the ruler integrates between the polygon's vertices. Another
device's photomask under the build's name fails P0, R0, R6b, P4, and P5, and M14 places it outside the envelope.
Tripling the flow puts the device outside the envelope (M14) without failing it. A lost feature in the DRC fails P9, a
corner sharper than the criterion fails P10, raising the target F by 0.01 fails P6, and an aspect-ratio limit below the
drawing fails P9-AR.

Finally, every shipped photomask was audited (`docs/VALIDATION.md`). {{AUDIT_SUMMARY}}

## Running the audit

```
python -m opencif audit runs/my_device                 # every check, its result, and its role; exit 0, 3, 2, or 4
python -m opencif audit runs/my_device --out audit.json
```

```python
from opencif.checks import auditor
a = auditor.audit_run('runs/my_device')
print(a['status'], [(c['id'], c['ok']) for c in a['checks']])
```

## The design code's other requirements

The auditor is one of the design code's sets of checks. The others that concern a design are its gates, which have
been ported as gates since rc3 (`build.gates` and `pipeline.prepare`), the mandatory reports of its `design()`, and the
operating validators of its acceptance run. Since 1.0.0, every build reports the last two by the design code's own
formulas, as the design code does, without making them pass conditions (`validation.design_code_reports`,
`opencif/checks/design_code.py`, and `tests/test_design_code.py`). `docs/CIF_V2_REQUIREMENTS.md` lists every gate,
validator, and acceptance criterion of `reference/CIF_v2.py` and its status in opencif 1.0.0.
