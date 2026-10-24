# opencif: the reference

This page is the full reference of the generator (RC8's README, kept current for 1.0.0), and `README.md` is the
starting page. opencif is a target-first generator for Controlled Incremental Filtration (CIF) devices. The user gives
the size cutoff D and the filtrate fraction F (two independent targets), the fabrication process, optionally a chip
area, and a choice of parts, and the generator returns a DXF photomask, three proof figures, and a validation record.
The record measures the drawn mask (every side segment, head unit, U-turn side and center arc, outlet lead, gap, and
center width) and reports, gap by gap, the cutoff that the 1-D network predicts for that drawing and the F that it
delivers.

One call builds a device from the two targets, as the following examples show.

```python
import opencif
r = opencif.generate(4.0, 0.5, out='/tmp/cif')    # D = 4 µm, F = 0.5; every other part defaulted (see docs/OPTIONS.md)
r = opencif.generate(2.5, 0.6, process='dinh2024', chip=(25, 75), layout='serpentine',
                     taper=('constant_shear', 70, 50), meander='pack', inlet='filter', outlet='teardrop', out='/tmp/cif')
print(r.verdict, r.exit_code, r.dxf, r.figures)   # r.record: the full validation record
print(r.ok, r.over_chip_area, r.over_aspect_ratio) # ok = exit 0; the chip area is reported, never part of ok
print(r.record['assessment']['status'])            # the verdict's parts as separate fields
```

`docs/OPTIONS.md` is the options sheet. It gives every part style, the code to write for it (YAML key and `generate()`
argument), its default, its origin, an example spec, and the record field that verifies it on the drawn mask. The same
builds run from the command line or from a YAML spec, as follows.

```
python3 -m opencif generate 2.5 0.6 --process dinh2024 --chip 25x75 --layout serpentine --taper constant_shear \
        --width-in 70 --width-min 50 --meander pack --inlet filter --outlet teardrop --out /tmp/cif
python3 -m opencif build examples/gifford_targets_meander_head_19um.yaml --out /tmp/cif   # default: runs/<name>/
python3 -m opencif check examples/platelet_serpentine_chip.yaml              # plan, gates and the outline's collision checks; nothing written (exit 0, 3 or 5)
python3 tools/make_arrays.py m07_serpentine_const_pack_filter_separate 3 shared 25 75   # 3 x a passing element, one inlet, 25 x 75 mm slide
#   (needs the element's runs/<name>/<name>.pkl from a build on this machine; it writes runs/<name>_x3/)
python3 -m opencif generate --help                                           # every argument, its units and choices
python3 -m opencif --list-processes                                          # the named processes and their sources
python3 -m opencif audit runs/<name>                                        # the physics audit again, from the photomask file (1.0.0)
python3 tests/test_kernels.py                                                # (1.0.0) the kernels = CIF_v2.py's, bit for bit
python3 tests/test_auditor.py                                                # (1.0.0) the auditor against CIF_v2's own; seeded faults
python3 tests/test_design_code.py                                            # (1.0.0) the design code's reports against its own functions
python3 tests/test_teaching.py                                               # (1.0.0) the teaching module = the generator's plan
python3 tests/test_units.py                                                  # unit checks (a few minutes)
python3 tests/test_presets_match_replicas.py                                 # presets = the replica package's designs (a round
                                                                             # trip to its factories, not to the papers' masks)
python3 tests/test_lattice.py                                                # the option lattice at the plan stage (10-20 min)
python3 tests/test_appendix_trace.py                                         # Appendix A/B results and glossary values pinned to the code
python3 tests/test_negative.py                                               # seeded defects: the request-side validators and refusals shown to catch them
python3 tests/test_bend_fd.py                                                # the curved-duct series against finite differences
python3 tests/test_options_sheet.py                                          # the options sheet covers every option; generate() reaches every part
python3 tests/test_figure_anchors.py                                         # presets and drawn masks against scale-bar measurements on the papers' figures (rc5;
                                                                             # exit 4 INCOMPLETE without the runs/ pickles; --source-anchors-only: the checks that need no build)
python3 tests/test_trace_anchors.py                                          # the generator's masks on the papers' micrographs, re-registered
                                                                             # from the rasters (OPENCIF_PAPER_RASTERS; provenance/rc7_traces/, rc6's in rc6_traces/;
                                                                             # RC8: coverage and censored nulls, and V13 on the Gifford Fig. 4 readings)
python3 tests/test_rc7.py                                                    # (rc7) the rc7 options and fixes, rc6 invariance of arc-turn heads (a few minutes)
python3 tests/test_composition.py                                            # (rc7 final) the composition sweep: a sample of layouts x parts; colliding layouts refused before writing
                                                                             # (~65 min; OPENCIF_TEST_CASES=0-10,x and 11-21 in parallel: ~35 min)
python3 tests/test_rc8.py                                                    # (RC8) input validation, the conflict pairs, the plateau pre-flight, the plan
                                                                             # notes, name_mask, the arrays' inter-copy rule (exit 4 without the runs/ pickles)
python3 tools/contact_sheet.py                                               # docs/CONTACT_SHEET.pdf: every layout, verdict and drawn parts
```

A build writes its outputs in `--out` (default `runs/<name>/`). It first removes that name's earlier files there, so in
an evaluation copy always give `--out` or `out=` to keep the shipped records. The outputs are the following.
- `<name>.dxf`, in which layer CHANNELS is the fluid and layer ISLANDS the posts.
- `<name>.pkl`, which holds the drawn polygon and the arrays from which the record was computed (used by
  `pipeline.reassess` and the tools).
- `<name>_audit.json`, the DRC audit of the written DXF.
- `<name>.json`, which holds the design, every check, the notes, the values taken from defaults, and the hash of the
  code that produced it.
- `<name>_summary.md`.
- Three figures:
  - `_layout.png`, the whole mask with the zoom windows marked.
  - `_zooms.png`, with the inlet, head, seam, first rows, a side-channel close-up, mid-array, U-turn, last rows, and
    outlet. The zooms are outlined, so thin channels stay visible.
  - `_physics.png`, with the widths, the drawn and design gap fractions, the realized cutoff against D over the whole
    device, at the entrance, and around the first U-turn, and F.

## How the two targets enter

* D sets the lamina width δ* = D/κ_g (cubic kernel pair, κ_g = 4.9138). Every gap therefore extracts
  f(i) = g(δ*/w_c(i)), so its cutoff is D whatever the center width. Nothing is computed by hand.
* F sets where the march stops, which is the number of gaps N (`tests/test_units.py`), with three consequences.
  - For a sealed start, changing D at fixed F changes N and leaves F unchanged.
  - Changing F changes N and leaves the per-gap cutoff unchanged.
  - A floor start stops on the delivered F.
* The closure is the exact gap-drop field rather than pressure equilibration. Every straight side-channel wall is drawn
  so that its integrated resistance is exactly the closure's, for which the midpoint widths are solved. The older
  linear rule left errors of up to 15% in f at a sealed start, and a midpoint-only check could not see them.
* Bends and junctions carry their own terms. Every arc carries the exact Stokes curved-duct correction
  (`schedule.bend_ratio`), less its two straight-arc junctions, and every meander-head unit ends in a T-junction counted
  with its own correction.
  - For the bend correction, the flow in a rectangular duct curved in its plane is purely azimuthal at Re → 0 and is
    solved in closed form by a Fourier-Bessel series. The correction is checked against its planar-annulus and
    Hele-Shaw limits and against an independent finite-difference solve.
  - The junction terms come from 2-D Stokes solves in the deep and shallow limits. A bend's junctions remove κ_j·w of
    the arc from the correction (κ_j 0.50 deep … 0.33 shallow). A head T-junction (a unit's drop turning into the next
    unit's stem) is shorter than the centerline convention by 0.28 w (deep) … 0.47 w (shallow) for sharp corners, plus
    a fillet term (r_f/w) × (−0.83 deep … −0.33 shallow) for the drawing's fillets, which is about a third of the
    term at the head's widths (28 µm units, 5 µm fillets: −0.29 w sharp, −0.14 w fillet). A depth weight taken from the exact curved-duct
    solution places a real channel between the limits (h/w = 5: 5% of the way to shallow; h/w = 1: 33%). These terms
    are fits rather than 3-D results, and their 3-D calibration is future work (`tools/cfd_kit.py` writes the windows).
  - The outlet stretch's end (rc5) carries a term of its own. Where the last side channel meets a lead of another width
    (a one-sided step into a wider lead, a centered step into a narrower one, or a flare), the flow's adjustment adds
    an excess of c·w of the narrower channel. 2-D solves give c in the deep (planar Stokes) and shallow (Hele-Shaw,
    checked against the conformal closed form) limits, by geometry and width ratio (and angle), and the same depth
    weight places it between them (`schedule.outlet_exp_c`; a right-angle step of ratio 2.4–2.8 is about 0.5 w, and a
    14.5° flare 0.08 w). In rc4 the term was left out, and on rc4's masks it was 0.34–3.49% of the side outlet's
    resistance, to which the last rows are very sensitive (LIMITATIONS 15a). The term applies to every outlet whose
    side channel changes width, so rc5 changes the outlet closure of every teardrop and ring design as well as Dinh's.
    Thus, 21 of the 46 design masks differ from rc4's at their outlet ends (the two Dinh presets also draw the measured
    flare), whereas the other 25, both arrays, and the five CFD windows are byte-identical
    (`provenance/RC4_REVIEWS_VERIFICATION_2026-10-11.md`). The term's band (the depth weight moved, c ±30%) is a
    sensitivity scenario rather than a calibrated interval (LIMITATIONS 8 says what it was checked against).
  - The terms apply in the design and in the measurement, to the U-turns (center and both sides), the head units'
    hairpins, arms, and T-junctions, the transition unit, the outlet leads, and the outlet stretch's end.
  - The terms matter because, without the bend correction, Gifford Fig. 4's turn rows move by up to 5.2% of D. A
    meander head's seam rows are the most sensitive in a headed device. Without the T-junction term, m09's seam rows
    move by up to 65% of D (54 rows outside ±1%), and a few tenths of a percent in the head's resistance move them by
    percent. The record's model band reports the rows that the junction terms' range can still move, with their
    share of F.
  - `physics.bend_model: none` restores the plain 1-D rule, for comparison only (evaluate mode; such an element cannot
    be arrayed).

## Input file (YAML; every key is validated, unknown keys or values raise an error)

| section | keys (units) | notes |
|---|---|---|
| `targets` | `cutoff_um`; `filtrate_fraction` F, or `filtrate_to_retentate_ratio` r = F/(1−F) (not the glossary's FR) | required |
| `process` | required: `depth_um` (h), `gap_um` (G), `min_feature_um`. Optional: `post_length_um` (along the flow, the glossary's g_p, "post width"), `post_width_um` (across the flow, the glossary's ℓ_g, "link length"), `min_radius_um`, `fillet_um`, `gap_entrance_c` (c_ent), `aspect_ratio_max` (AR_max, default 10), and, since rc7, `wall_along_flow_um` / `wall_cross_flow_um`, the glossary's names for the two post dimensions (aliases, which must agree when both are given) | The floor is required, because every check runs against it. AR_max gives CIF_v2's w_min = max(floor, depth/AR_max) for its width gates and for the aspect-ratio check on the drawn mask (a pass condition in strict mode, reported in the other modes). |
| `chip` | `width_mm`, `length_mm`, `margin_um` | optional; enables leg-count fitting (and head widening) and the footprint report (a design over its chip is drawn and reported and is never refused or failed, App. B C8 CQ5) |
| `centre_channel` | `taper`: `constant` (Gifford) or `constant_shear` (Xia; Dinh's form assumed); `width_in_um`; `width_min_um` (constant_shear only); since rc7, `shear_formula` (constant_shear only), which is `exact` (default, C(α)) or `parallel_plates` (6 Q/(a b²), reproduce mode only, for Xia's figures) | |
| `layout` | `style`: `linear`, `serpentine` (odd legs), or `paperclip` (nested, as in Gifford Fig. 4 and Xia). `legs`: an integer or `auto`; since rc7, `order` (paperclip), which is `nested` (default) or `path` (legs side by side in path order) | `auto` gives the fewest legs that fit the chip and draw without collisions (one polygon with the planned holes, and outlet ports separated from each other and from the inlet; a two-leg paperclip ends beside its inlet). If no count fits, a packed head with arms is widened (its arms fold into more lanes) and the fit is repeated, whereas hairpin units, Dinh's arms, and the narrow head do not shorten with width. |
| `meander` | `enabled`; `style`, which is `pack` (the default, Dinh's arms topology with as many arms, folded into lanes, as the head width holds, then hairpins), `hairpins` (Xia 2016), `arms` (Dinh 2024; `arms: n`, single lanes), or `narrow` (opt-in, with process-floor units first, a smaller head, and three times the width sensitivity); `channel_um`, `wall_um`, `narrow_um`, `head_width_um`; `seal`, which is `require` (gap 1 sealed by the asked head, or the build is refused with the numbers, which are the arms needed and the room for them and the head width that would seal; no other head is substituted) / `auto` (as require in strict mode; otherwise an unsealed gap 1 is declared and reported) / `allow` (the asked head unsealed and gap 1 declared; not accepted in strict mode); `head_width_um` (by default the narrowest width that seals, searched from 0.3 to 1 times the chip width less the margins); `transition`, which is `exact` (`paper` is reproduction-only); `rectilinear_rows` (≥ 1, default 2), the wall-block rows between the head and the stadium posts; and, since rc7, `turns`, which is `arc` (default) or `square` (right-angle quarter turns, `corner_radius_um`, `inner_corner_radius_um`, drawn at no less than the fillet), and `junction_radius_um` (the wall corners where each unit leaves the center channel; drawing only) | a meander replaces the narrow start of the side channels with equal-resistance segments (Xia 2016) |
| `start` | `style`: `sealed` (default) or `floor` (Gifford's rule, by which the first side channel is drawn at `side_width_um`) | no-meander devices only. `floor` over-extracts at the entrance. The design network predicts those rows, and the drawing is measured against that prediction. It stops on the delivered F. |
| `inlet` | `style`: `plain`, `filter`, or `dinh_filter` (Dinh's, in absolute coordinates); port radius, stem | `filter` draws an aggregate filter, a bulb of half-width up to 150 µm with staggered 30 µm posts (60 µm pitch), between the arms or above the head, and it narrows to the trunk 50 µm before anything else is drawn. The filter lies upstream of gap 1, so it does not enter the network. |
| `outlet` | `style`: `separate` (Gifford Fig. 4), `teardrop` (Dinh), or `ring` (after Xia 2016 Fig. 3a); port radii; `stretch_um`; `flare_deg`, `flare_start_um` (rc5); `tubing` (I.D.s, length); `ring` geometry | The ring is parametric (scaled on the ports) unless given, and it is available in linear and paperclip layouts only. `flare_deg` (1–45°, per side; teardrop and separate only; refused with `stretch_um` and with `physics.fail_safe_margin`) draws the side channel straight to `flare_start_um` (default 67, the flare's virtual corner) from the last gap's center, or farther where the last leg ends beside a loop, and then flares its outer wall to the lead width, as in Dinh 2024 Fig. 2B. A lead no wider than the side channel draws no flare, and when `flare_deg` is unset the side channel steps into its lead at a right angle. The record's `design.outlet` gives the drawn stretch (`L_out_um`) and its end term (`stretch_end`). An outlet that cannot close the design (a side lead that would need less resistance than its straight run has, or one too wide for the ring's turn) is refused by name before drawing, in every mode, because it is planned on every leg count that the layout would try (rc4; exit 5, `REFUSED (outlet)`). `check` plans it at the head width that the spec gives, and the build plans again a head that the fitter widens later (LIMITATIONS 18). |
| `operation` | `flow_uL_min` | used by the Reynolds and Stokes gates |
| `fluid` (rc7) | `density_kg_m3`, `viscosity_Pa_s`, `rheology: newtonian` | default water (ρ 1,000 kg/m³, μ 1 mPa·s); used by the Reynolds and Stokes gates only, because the schedule's resistance ratios do not depend on μ, and the record's absolute resistances (Pa s/m³) stay water values |
| `gates` | `mode`: `strict`, `evaluate`, or `reproduce` | `strict` refuses a design outside the validity envelope, whose limits are the κ plateau and foot; channel Re ≤ 10 and the Stokes gates Re_gap and Re_δ ≤ 10⁻², as App. B's thesis states of the design code (CIF_v2 M14); and gap dominance M13 < 0.05, thin lamina E11 on every row, and taper rate M16 (CIF_v2 M13, E11, M16, which are not App. B gates). It also refuses a design that fails CIF_v2's width and γ_GS gates, checked before drawing (γ_GS; G and GS before planning, and w_c and every straight side channel w_s ≥ w_min once planned, by spec v5 G7, where G11 would fold a narrower side channel; rc4). On the drawn mask, strict mode fails a CIF channel above the aspect-ratio limit (exit 2, a failure rather than a refusal), whereas `evaluate` and `reproduce` report it. `evaluate` builds and flags. `reproduce` builds and allows the papers' declared departures. |
| `checks` | `cutoff_tolerance` (default 0.01), `F_tolerance` (default 1e-3), `drc`; since rc7, `passage_report` (default on, a report of each channel class against D, 2 D, and the largest realized cutoff) | the tolerances are at most 0.05 and 0.01 in every mode |
| `physics` | `fail_safe_margin`, `trim_F`, `wall_rule`, `bend_model` | The defaults are off, off, `exact`, and `stokes`. `bend_model: none` is refused in strict mode. |
| `array` | `count`, `pitch_um`, `inlet` (`shared` / `separate`), `bus_width_um` | a parallel array of the element (see Composition) |
| `reproduction` | engine-level overrides for drawing a published device as published | `reproduce` mode only. Every override is listed in the record. |

## What every build checks

* The fabrication checks run on the drawn polygon.
  - The DRC runs at the floor, tested at floor − 0.1 µm. It must run, or the verdict is INCOMPLETE.
  - The topology check requires every post, unit, and loop to be present.
  - Side-channel continuity requires one strip per leg and side, with the narrowest straight or U-turn section checked
    against the floor.
  - The outlet ports are separated.
  - No parts intersect. Every part (the inlet with the head, each leg, each U-turn, and the outlet) lies at least the
    floor from every part that it is not joined to, measured on the drawn outlines and named in the record; the outlet
    ports are separated; and every outlet port is at least 2 mm (edge to edge) from the inlet port. The leg fitter
    skips any leg count that breaks these conditions (when every count does, it hands the build the first whose outlet
    alone meets the inlet end). Where only the outlet meets the inlet end, the build lengthens the outlet stretch and
    re-plans (recorded as `design.outlet_stretch_extended_um`, and not done for a paperclip's ring outlet). Since rc7,
    the conditions are checked on the outline in memory before the photomask is written. Parts that collide are
    REFUSED (exit 5) with the two parts named, the remedies (for a paperclip whose spacing or leg axes the spec gives,
    the leg spacing that the generator's own rule gives), and, when figures are drawn, a figure
    (`<name>_refused_layout.png`), and no photomask exists. A spacing or leg axes that the spec gives are used as
    given, so this is where a spacing too tight for the head is caught. `opencif check` and
    `generate(check_only=True)` draw the same outline in memory and refuse the same way.
  - The drawn components are those requested. The check compares with the spec, rather than with the plan, the head
    style (hairpins with no arms and one width, Dinh's arms in the asked number and in single lanes, a packing head with
    full-width units, and a narrow head with units at the narrow width), gap 1 sealed when `seal: require`, the filter
    posts counted on the drawing for a filter inlet (none for a plain one), the outlet style and port count, and the
    leg count. No head is ever substituted. A head that cannot seal gap 1 is refused with the numbers (the arms needed
    and the room for them, and the width that would seal).
  - The head envelope.
  - The corner radius.
  - The chip area (when a chip is given) is reported and is never a pass condition. A design over its chip is drawn,
    and its verdict ends "OVER THE CHIP AREA ... reported, not refused or failed" (App. B's thesis states that the
    design tool "reports, without refusing, a design over its length budget"; C8 CQ5; the 2-D chip footprint
    generalizes that 1-D length budget). `assessment.over_chip_area` and `Result.over_chip_area` say so apart from the
    exit code.
  - Spec v5 G9's own length budget, which compares each leg's path length with the chip's leg length, is reported
    (`validation.length_budget`, `status.length_budget_G9`; rc4).
  - The aspect ratio of the narrowest CIF channel, measured on the drawing (gap openings, center channel, side
    channels, U-turn sides, head units, and outlet leads; rc5), is checked against `process.aspect_ratio_max`. It is a
    pass condition in strict mode (spec v5 G7; rc4) and is reported in evaluate and reproduce modes, where the verdict
    ends "ABOVE THE PROCESS ASPECT-RATIO LIMIT ... reported, not refused or failed in ... mode"
    (`assessment.over_aspect_ratio`; LIMITATIONS 7).
  - `assessment.status` gives the verdict's parts as separate fields (the hard checks, the envelope, the chip area, the
    aspect ratio, the G9 length budget, the drawn-mask assessment that was run, the physics audit of the photomask file
    (`status.physics_audit`, which is PASS, PASS_OUTSIDE_ENVELOPE, FAIL, or INCOMPLETE), the statement that nothing
    was fabricated or measured, and, for a reproduction, what the label covers).
* The physics audit (1.0.0) runs on the photomask file read back. It is the design code's physics auditor (CIF_v2's P1
  to P10), ported, with P0 (the file is the drawing), R0 (the auditor's own ruler re-measures every ordinary straight
  row from the file, and the generator's measurement must agree within 10⁻³), R6b (every straight-row post is the
  process's post), and the design code's gate M14 on the drawn mask. Its hard checks are pass conditions in every mode
  (exit 2), M14 is part of the validity envelope (exit 3), and P7 and P8 (CIF_v2's stricter inlet limit) are reported
  (`validation.physics_audit`; `docs/PHYSICS_AUDIT.md`). `python -m opencif audit RUN_DIR` repeats the audit from
  scratch.
* The design code's reports (1.0.0) are the mandatory reports of CIF_v2's `design()` and the operating validators of
  its acceptance run, computed by its own formulas and recorded, and they are never a pass condition
  (`validation.design_code_reports`).
* Q1 asks whether the mask draws its design. The 3-rail network is solved on resistances integrated over the drawing
  and compared with the intended network (the closure, plus any declared deliberate deviation), within the following
  limits.
  - At every gap, the drawn and intended cutoffs differ by at most 0.25% of the larger of the intended cutoff and D (a
    deliberate departure far above D is judged on its own size).
  - Where a row has no cutoff (reverse flow, f ≥ ½), f is within 0.5% of the larger of |f| and f*.
  - F is within 1e-4.

  The intended network of the papers' rectilinear transition is the profile laid down by the drawing code, and it is
  never copied from the measurement. It is integrated with the same quadrature as the measurement.
* Q2 asks whether the design realizes the targets.
  - Every gap's realized cutoff is within D ± tol on both sides, and rows with reverse flow or f ≥ ½ count as failures.
  - F is within F_tolerance of the target.
  - In `reproduce` mode, a row is exempt only when a declared option is shown to move it. Each row that the intended
    network puts outside D ± tol is attributed before measuring, in one of two ways.
    - It is attributed to the papers' transition or the pull-in when the intended network without that option brings
      it back inside.
    - It is attributed to the seed (gap 1 only), a depicted head (its own rows only), or a floor start (its leading
      run of entrance rows).

    A row that no option explains is not exempt, whatever its size. Departures are reported per option.
  - In `evaluate` mode, an explicitly allowed seed exempts gap 1 only, and the verdict is "PASS with an allowed
    departure".
  - Nothing is exempt otherwise.
* The validity envelope (App. B thesis scope; CIF_v2 M13/M14) has the following parts, and strict mode refuses a design
  outside it before drawing.
  - The κ plateau u = D/w_c ≥ 0.035 (the glossary's d̃, an external input) and its foot 0.012. Below the plateau, the
    record gives the cutoff band by App. B C7's method (the frozen Hill fit, plateau-anchored; App. B C7 CQ3).
  - A channel Reynolds number of at most 10 at the stated flow.
  - The Stokes gates, Re_gap ≤ 10⁻² and Re_δ ≤ 10⁻² at every row (CIF_v2 M14, its row-1 term included).
  - Gap dominance M13 < 0.05 (spec v5 §19).
  - The thin lamina, δ(i)/h ≤ 0.10 at every planned row (spec v5 E11 "∀i", since rc4; CIF_v2 checks the cutoff rows,
    whose lamina is δ*), and the taper rate, |w_c(i) − w_c(i−1)|/(G + g_p) ≤ 10⁻² at every row (CIF_v2 M16), with
    CIF_v2's tolerances.
  - The fluid is water unless the spec's `fluid` section says otherwise (rc7; ρ 1,000 kg/m³, μ 1 mPa·s). The
    quantities reported rather than gated are the largest f_gap against κ_g's calibration cap 5.76 × 10⁻⁴ and spec v5's
    f_soft 5 × 10⁻⁴; L_e/Λ (the glossary says "do not gate on it"); the quantities of CIF_v2's physics auditor P8
    (inlet Re against 1.0, and 0.05 Re D_h against Λ/2), at the design flow and at CIF_v2's default 5 µL/min
    (LIMITATIONS 6 explains that a P8 gate would refuse every design at its own flow); and a Reynolds bound for every
    channel, 2ρQ/(μ(h + F_min − 0.1 µm)), with Dean's leading term at it where that bound is ≤ 10 (LIMITATIONS 11).
* The 1-D model band is reported. For it, the same drawn network is re-solved in four ways.
  - It is re-solved with the junction terms (bend, head T-junction, and, since rc5, the outlet stretch's end, together)
    at both ends of their range, with the depth weight moved by ±max(weight, 0.1), the T-junction term widened by
    ±(a quarter of its fillet part + 0.02 w), and the stretch-end term by ±30%.
  - For information, it is re-solved with the T-junction term ignored, with all junction terms ignored, and with no
    bend correction.
  - It is re-solved with both side outlets' resistance ±1% (the outlet balance on which the last rows rest,
    LIMITATIONS 15a; rc4; reported in the record rather than in the verdict).
  - It is re-solved with the off-chip tubing varied (length ×0.5 and ×2; center I.D. ±0.025 mm).

  Changes are counted on the rows that are inside D ± tol on the drawn mask. When any of them leaves D ± tol, the
  verdict says so in brackets, with the rows. Dean (inertial) terms are not modeled, because Dean's leading term is
  about 1.5 × 10⁻⁴ of a bend's resistance at Re = 10 and a/R = 0.2, and 3.7 × 10⁻³ at a/R = 1, the most that any
  drawable arc reaches. A record gives that screening estimate (a circular pipe's term rather than a validated bound)
  for every bend of its device where its any-channel Reynolds bound is ≤ 10 (LIMITATIONS 11).
* Scenarios of the design network are reported. They are a head-unit junction length ±w/2, a uniform width bias
  ±0.5 µm, the head mouths' R_gap ±20%, and one seeded draw of 10 nm per-segment side-width noise (a single
  realization; there is no depth scenario). The sizes are illustrative assumptions rather than measured process
  tolerances. The scenarios perturb the design network rather than the drawing, and they are never a pass condition
  (see Limits for their size).

The exit codes are as follows.

| code | meaning |
|---|---|
| 0 | PASS or REPRODUCED |
| 2 | FAIL (a hard check, the physics audit's included); the photomask is kept for diagnosis as `<name>_FAILED.dxf` and never under the deliverable name `<name>.dxf` (rc7) |
| 3 | outside the validity envelope in evaluate or reproduce mode, and drawn and measured (build) or planned and outlined in memory (`check`, `generate --check`; nothing written) |
| 4 | INCOMPLETE (for example, no DRC) |
| 5 | REFUSED, with no mask written. Before drawing, a design is refused in strict mode outside the envelope, or when it fails a gate checked before drawing (spec v5 §14.4, §18, G7; CIF_v2 `design()` and its validator G7), which covers a meander head's GS/G outside γ_GS [1.4, 1.5] and the gap, a drawn GS, the center channel, or any straight side channel below w_min = max(floor, depth/AR_max). In every mode, a design is refused for parts of the drawn outline that collide (since rc7, found on the outline in memory before the photomask is written, and the record names the parts and the remedies, with a figure), for a head that cannot seal gap 1 under `seal: require`, or for an outlet that cannot close the design (rc4; found before drawing, or while fitting or drawing). The cutoff band is still reported. |
| 6 | SPEC ERROR (command line), for an invalid or unplannable spec (`generate()` raises `opencif.spec.SpecError`). Since rc4, its message also names every other problem found before drawing. Since RC8, a negative, zero, infinite, or true/false length, depth, flow, radius, tolerance, or count (nine keys may be 0), an unknown process, and contradictory head widths are spec errors with the field named. |
| 1 | the generator could not draw the request (an exception that the checks above did not anticipate), and no record is written |

## Composition

* In series, within one device, the parts run head → transition → ladder legs → U-turns → outlet. Every part is
  resistance-exact and checked gap by gap.
* In parallel, `opencif.compose.array` (driven by `tools/make_arrays.py` or by the `array` block of a spec) stacks n
  copies of a validated element.
  - Only a clean PASS element can be arrayed, which means exit code 0, inside the envelope, and no exempted rows.
  - Each copy keeps its own outlet ports and the tubing with which it was validated. The 1-D network is linear (its
    resistances do not depend on the flow), so every copy then realizes the element's D and F whatever the inlet
    manifold does.
  - The inlet manifold only splits the flow. The split is computed, and the Reynolds gate is applied at the busiest
    element.
  - A shared inlet bus is drawn only when every inlet port can be reached in one layer, as in a linear or serpentine
    element. Otherwise, each element keeps its own inlet.
* Merged outlets on the chip and multi-stage cascades are not implemented. Merging outlets changes each element's
  outlet closure.

## Presets (reproductions) and examples (new designs)

* `docs/PRESET_PROVENANCE.md` gives the source of every preset value, which is printed, measured, derived, assumed, or
  fitted.
* `docs/RESULTS.md` gives each preset's and example's current verdict, as produced by the code hash in its record.

## Limits

`docs/LIMITATIONS.md` lists every limitation and the future work that it implies. `docs/ALIGNMENT.md` maps every
statement of Appendices A and B and the glossary to the code and the test that pins it. `docs/REPRODUCTIONS.md` puts
each paper's printed quantities beside the reproduced ones. `provenance/` holds the scripts behind the values that the
package takes as inputs (the junction fits, the figure measurements, and the Fig. 5 reading). `docs/EVALUATION.md` is
the guide for a third party evaluating the code base.

* The scope is rigid, approximately spherical particles in a Newtonian, nearly constant-viscosity fluid, with laminar
  channel flow (Re ≤ 10, gated) and creeping flow at the gap and lamina scales that select particles (the glossary
  states "The framework is therefore not globally Stokes", and App. B's thesis "rigid spheres in Stokes flow").
  Deformability, cell crowding, and whole blood are out of scope, and platelets are discoid, so platelet-specific
  validity is not established. D is the rigid-sphere cutoff.
* The cutoff is a model prediction. The realized cutoff is the 1-D network model evaluated on the drawn geometry rather
  than a measured particle cutoff. κ_g and the plateau come from Gifford 2014's calibration. 3-D flow at gap mouths and
  at the outlet bifurcation is outside the model.
* D ± 1% at every gap is a statement about the nominal geometry. Per-gap cutoffs are very sensitive to sub-micron width
  errors. On the shipped designs, an illustrative uniform ±0.5 µm width bias moves 17 to 223 gaps of a device without
  a head (up to 21% of N; cutoffs 0.84–1.16 D, worst at a sealed start's first rows) and 82 to 404 gaps of a headed
  device (up to 23%; declared rows, outside at nominal, counted too; worst at the seam rows; counts over both signs).
  One seeded draw of 10 nm per-segment noise moves 0 to 2,005 gaps. Every record carries these numbers
  (`validation.robustness_design_network`), and RESULTS.md lists them per design (LIMITATIONS 13–14).
* The bend correction is the fully developed Stokes solution, whereas the bend, head T-junction, and (rc5) outlet
  stretch-end terms are fits to 2-D solves in two limits, placed between them by a depth proxy. The head's seam rows,
  the rows next to a tight U-turn, and the last rows (which rest on the outlet balance, LIMITATIONS 15a) are the most
  sensitive in the device (see the model band).
* Some parts are taken from the design rather than measured.
  - The radii of the head, U-turn, and outlet arcs (their widths are measured).
  - The head rows' center-path lengths.
  - Any off-chip tubing.

  Since 1.0.0, the straight rows' lengths and the post length across the flow (l_g) are measured on the file by the
  physics audit's own ruler (R0) and its post check (R6b) (LIMITATIONS 12).
* The tapers are constant and constant-shear only, and there is no linear, canonical Ψ (App. A R13(b)), or user-defined
  profile (LIMITATIONS 16).
* Composition is described above.
