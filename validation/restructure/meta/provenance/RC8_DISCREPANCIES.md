# opencif RC8: the discrepancy register (rc7's register, with the RC8 entries)

Started 2026-10-18, while building the four device packages on the frozen rc6 code (stamp `4bbdaabef1225a24`).

**RC8 (2026-10-21; code `593d977fc83ec217`).** Two third-party reviews of rc7 (2026-10-20), their adjudication
(`RC7_REVIEWS_ADJUDICATION.md`), a third review of the RC8 candidate (`RC8_CANDIDATE_REVIEW.md`) and its adjudication
(`RC8_FINAL_REPORT.md`) changed these entries; everything else stands as rc7 left it:

- **GIF4-1, INV-2 (decision changed).** The Gifford Fig. 4 reproduction is the nested reading (D1), under the same name
  `gifford2014_fig4_aligned` (D6): the inset's taper direction rejects reading B, rc7's device, which is kept as the
  variant `gifford2014_fig4_pathorder`.  Leg 2's side-channel widths are not reproduced (open).
- **XIA-9 (new).** The trace metric's coverage (A7): Xia Fig. 3d (both masks) and 3g (the depicted mask) are recorded,
  not claimed (D5).
- **BUILD-2 (closed further).** The contradictory-input pairs, the array's stale photomask, the extension's tries and the
  arrays' inter-copy port rule (K7, K8, K10, K11) are fixed in the code.
- **XIA-1 (note).** Xia's sentence on the taper's method is quoted in `docs/ALIGNMENT.md` S1 as well.

**The plan (A. H. Routt, 2026-10-18 14:09).**

1. Work through the four device packages on rc6, one at a time.
2. Note every discrepancy between a paper's figures and the reproduction, with its possible solutions.
3. Then make rc7 the final frozen codebase: verify, run independent reviews until they converge, freeze.
4. Re-run the four packages on rc7. Each package is one script, so this is one command per device.

**Rules.**

- **Evidence.** The only evidence is the papers' published figures and text. No mask files or other device data are available, so whatever the figures do not show is noted, not drawn.
- **Kinds of entry.** Each entry says what the figure shows, what the mask draws, the size of the difference, the cause if known, the possible solutions, and how each solution stays within:
  - the frozen foundation (Appendices A and B, glossary v114);
  - the design specification (spec v5);
  - the rc6 code.
- **Status.**
  - **rc7**: a code or record change planned for rc7.
  - **doc**: a documentation change, made at rc7.
  - **note**: recorded and disclosed, not changed. Either the figures cannot settle it, or it is a property of the paper.
  - **done**: fixed in the packages, with no rc7 change.
  - **open**: a decision is pending.
  - **RC8**: changed at RC8 (2026-10-21): a decision, a code or record change, or a document change.

Codes for values follow PRESET_PROVENANCE: **P** printed, **M** measured on a figure, **D** derived, **A** assumed, **G** the generator's choice.

**The aim (A. H. Routt, 2026-10-18 17:25–17:28): full alignment.** A reproduction aligns when the paper's figures, the frozen physics (App. A/B, glossary, the code) and the printed values all agree on it: every visible feature lies on the drawing, every width comes from the physics, f_gap and the final width (so F) are the printed ones. Where the paper contradicts itself (a caption against its text, as at several places in Gifford 2014), the reproduction follows the value the figures support, and the dissertation states the conflict. A start width, gap or other value Gifford did not print is measured on the figures (**M**) and labelled so; what a reproduced device then does (for Gifford's devices: the first rows filtering coarser than D) is reported as a result.

---

## 1. Device 1: Dinh et al. 2024, V1, Fig. 2B (rc6 package: `dinh2024_v1_depicted`; rc7 package: `dinh2024_v1_depicted` with square head turns, and `dinh2024_v1_figure_aligned`)

### DINH-1. The head's unit turns: arcs where the device has square corners (rc7)

**What the figure shows.** In Fig. 2B's first panel, each returning meander channel turns into the next branch at a square outer corner, with a small fabrication rounding. Measured on the panel at its rc6 registration (1.0660 µm/px), by the outline's distance from the apex of the two fitted straight walls (`provenance/rc7_scripts/measure_p1.py`):

- left: 6.0 µm;
- right: 8.6 µm.

The ends of the walls between channels (the hairpins' inner sides) are rounded:

- the best corner radius is 9–11 µm;
- a full semicircle (radius 12.3 µm, half the 24.5 µm wall read on the image) fits as well on the left (rms 1.60 against 1.58 px) and 28 % worse on the right (1.15 against 0.90 px).

**What the mask draws.** `engine/geometry.py`, `unit_path`: every turn is a constant-width arc of centreline radius s/2 = (28 + 27)/2 = 27.5 µm. That puts the outer wall at radius 41.5 µm and the inner wall at 13.5 µm. The same rule applies to the arm lanes.

**Size.** At a corner with rounding r, the arc lies (√2 − 1)(41.5 − r) from the corner along its bisector:

- 14.7 µm for r = 6.0;
- 13.6 µm for r = 8.6.

Near the apex that is beyond the metric's ± 8 px (8.5 µm) window. So those stretches drop out of panel 1's median (1.29 µm) and p90 (2.70 µm), and the metric found an outline on only 53 % of panel 1's classified samples (a lower bound). The wall ends' semicircles are within about 2 µm of the image.

**Cause.** An assumed shape (**A**), not a reading. The source of the head lengths, Fig. 2A's inset, draws every corner square. So the depicted lengths were measured on square-cornered channels and then drawn on arcs.

**Physics today.** The head's resistance is computed on the arcs drawn, with rc6's exact curved-duct correction (`physics.bend_model: stokes`). The model is therefore consistent with the mask, but not with the device's corners.

**Possible solutions.**

1. **(Preferred) A turn-shape option.**
   - **The option.** `meander.turns: arc | square` (default `arc`, today's drawing, so every existing photomask stays byte-identical) and `meander.corner_radius_um` (outer corner rounding; inner corners at least R_min).
   - **The Dinh preset.** `dinh2024_v1_depicted` would take `square` with r = 7 µm (M, the mean of 6.0 and 8.6; its band ± 1.3 µm), and the wall ends would keep their semicircles.
   - **Physics.** A right-angle turn's resistance is the straight centreline length plus a bend-equivalent length L_eq per turn:
     - spec v5 § 14.6: R_s^eff = r(GS)[L̂ + n_turns L_eq];
     - the glossary (v114, External Inputs table) lists L_eq, "the bend-equivalent length", as **external, supplied by this work**, and R_min as the minimum inner-corner radius;
     - Appendices A and B say nothing about meander shape.

     So no frozen text changes; rc7 supplies L_eq.
   - **Computing L_eq.** The same machinery rc6 already uses for its junction terms (`provenance/junction_fits/`): 2-D solves in the deep limit (in-plane Stokes) and in the shallow limit (Hele-Shaw), with rc6's depth weight between them, and a model band. The corner has its inner and outer fillets at the measured radii; rc6's fits already include a plain 90° corner (−0.24 w deep, −0.44 w shallow).
   - **Validating the solver.** The solver is first checked against rc6's exact curved-duct solution for an arc turn. If it does not reproduce rc6's bend correction to within its band, it is not used.
   - **Recording.** The record carries n_turns and L_eq per unit, so that spec v5 G12 ("n_turns counted from the drawn polyline") holds.
   - **Effect.** Each square quarter-turn's centreline is longer than the arc's by (2 − π/2)·(s/2) ≈ 11.8 µm, before L_eq. With the depicted head's lengths fixed (Fig. 2A), the head units' resistances change by a few per cent. The head's per-row fractions (declared rows 2–18, already outside D ± 1 % by design, "head drawn as depicted") move. The linear section is unchanged. The DXF changes: a new rc7 photomask for this preset.
2. **Keep arcs, declare only.** This is rc6's state, plus the disclosure in LIMITATIONS 24 (adjudication edit A6). It costs nothing, but leaves the panel-1 overlay visibly wrong at every turn.
3. **Shrink the arc radius** to the measured corner (s/2 → about 7 µm). This is impossible at constant width: a 28 µm channel cannot turn about a 7 µm centreline radius. Rejected.

**Status.** done in rc7: `meander.turns: square`, outer 7 / inner 10 µm (both rc7 Dinh presets that draw the depicted head), c_corner per corner (`model_scenarios.bend_model.square_corner`); n_turns is not recorded per unit: G12 holds by construction (LIMITATIONS 8).  The plan above ("r = 7 µm … the wall ends would keep their semicircles"; "Recording … n_turns and L_eq per unit") is not what shipped: the inner corners are 10 µm.

### DINH-2. A closed outline on each side, upstream of the first branch, that the mask does not draw (note)

**What the figure shows.** In Fig. 2B's first panel, above the first branch on both sides, there is a closed, rounded-rectangular outline. Measured at the registration:

- **Along the flow:** it spans x = −122 to −39 µm in mask coordinates (83.2 µm on the left, 83.5 µm on the right).
- **Across:** its end nearer the axis lies 105 µm from it (104.9 and 105.2 µm).
- **Corners:** rounded 5.7–8.2 µm.
- **Extent:** it runs out of the frame, more than 170 µm from the axis.

83 µm equals 28 + 27 + 28 µm (two meander channels and a wall), but no wall shows inside it within the frame.

**What the mask draws.** Solid wall (no feature) between the inlet funnel and the first branch.

**Why it cannot be drawn.**

- Fig. 2A's inset, the only drawing of the head, shows nothing there: its top line is the first branch itself.
- The outline leaves the frame, and no published figure shows its course.
- The panel cannot tell whether it is fluid (a channel, or a dead-end pocket) or a structure of the PDMS. The relief image shows an outline, not which side is fluid.

So drawing it would be a guess.

**What it would change, if it were known.**

- A dead end (fluid that leads nowhere) carries no flow: no change to the network or any cutoff.
- A part of the first arms' path would change the first units' side resistance R_s(1..3), and so the head's first fractions (declared rows).
- A PDMS structure would change nothing.

**Status.** note. A. H. Routt, 14:09: only the papers' figures are available, so this is noted, not reproduced. Disclosed in the Dinh package's README.

### DINH-3. Device length: 148.2 mm drawn (rc7, square turns; rc6 148.3), "∼3.2 mm × 75.8 mm" printed (note)

This is a conditional conflict under the model's premises, not a conflict among the printed numbers. It is declared in rc6: the preset's `paper.length_bound`, PRESET_PROVENANCE, REPRODUCTIONS, DISSERTATION_NOTES. The micrographs register locally, and their places along the device are the mask's.

No code change. The dissertation text carries the paragraph (adjudication 12-E).

### DINH-4. Transition blocks read about 2.2 µm short at each end along the flow (note)

Fig. 2B's second panel: V3 block_along −2.10 / −2.47 / −2.23 µm (centre / left / right). LIMITATIONS 24 (a). A relief image's outline bias and blocks drawn 2 (s − GS) long cannot be separated on one panel.

### DINH-5. Outlet: the centre channel's widening over the last ~37 µm is not drawn (note)

Fig. 2B's fourth panel. The lead width after the flare (105.75 µm) is the generator's closure (**A**): it is not visible in the panel.

### DINH-6. Fig. 3C is not traced (note)

At their printed 100 µm bars the post pitch reads 51.2–51.4 µm, against 45.1 µm on all three post-bearing panels of Fig. 2B at their 50 µm bars. So one of the two figures' bars is not the stated length. Fig. 3C's images ship in `inputs/` but are not used.  Resolved: INV-5 (Fig. 3C's 100 µm bars read 12–14 % short).

### DINH-P1. The package's input table covered 26 of the preset's 36 governing values (done)

The ten values it missed:

- `meander.enabled`, `meander.seal`, `meander.transition`, `meander.seam_at`;
- the three port radii and the stem half-width;
- `reproduction.margin`, `reproduction.ws_rect`.

Fixed at 14:13 on 2026-10-18: `devicekit.check_preset` now requires the table to hold every value of the preset (the free text under `name`, `description` and `paper` aside), with a code from P/M/D/A/G. Package only; no rc7 change.

---

### DINH-7. The head's corners on Fig. 2B: all rounded; the trunk-junction corners about twice the drawn fillet (rc7, measured)

**The question** (A. H. Routt, 2026-10-19 00:23): "all inner corners that impact particle flow are necessarily still rounded".

**Measured** (`review7/dinh1/scripts/corner_radii.py` on the rc7 square-turn mask's panel-1 registration; the rounding from the apex's distance to the outline along the bisector, d = (√2 − 1) r; `review7/dinh1/trace/corner_radii_p1.json`):

- every corner of the panel is rounded on the device;
- the wall corners where the head units' stems meet the trunk (the head's "gaps": the corners particles pass at extraction), 8 corners: r ≈ 11 µm (median 11.3; 9.0–14 µm, two readings of 24–33 µm beside the panel's arrows and label left out), against the 5 µm head fillet drawn;
- a quarter turn's inner corner (a wall's end): 7.3 µm against the 10 µm drawn (one corner); the outer (fluid) corners: 6.4 µm median of three noisy readings (the dedicated measurement of DINH-1: 6.0 / 8.6 µm), against 7 µm drawn.

**What it moves.** These corners are the head gaps' entrance from the trunk (y = w_c0/2 = 37.5 µm), not the T-junction where the gap meets its unit (y = w_c0/2 + l_g = 91.5 µm; its wall corners read 7.6 µm and one noisy 18–33 µm against the 5 µm fillet). They enter the model only through the gap entrance coefficient c_ent (external, 1.0 at every gap; the design is "insensitive while M13 passes", spec v5 S:756).  Their particle-scale effect is outside the frozen cutoff model: κ_app is calibrated on Gifford's array (G:741–743), whose post shape the foundation does not state (G:760–762), and the head rows lie below the plateau anyway (B:4146–4149).  The T-junction corners, which c_T's fillet part does model, read 7.6 µm against the 5 µm fillet: c_T would move by about −0.075 w (−2.1 µm per unit, 1.3 times its band); disclosed (LIMITATIONS 27), not drawn (one clear reading). (Corrected after the code review: the first note applied the T-junction's fillet fit, c_T −0.433 w at 5 µm to −0.531 w at 8.4 µm, to these corners; that does not apply.)

**Robustness.** The reading depends on the registration: on the arc-free square-turn trace the six corners away from the panel's arrows read 9.0–13.8 µm (median 11.4); on the trace of the mask that draws them at 11 µm, 6.9–14.5 µm (median 14.0), with two to three readings of 22–33 µm beside the arrows and label. So: larger than the 5 µm fillet, about 7–15 µm; 11 µm is a mid reading, not a precise one.  Rerun on the final rc7 trace (code 4a624626f570918c, the shipped mask; `provenance/rc7_scripts/dinh_corners/corner_radii_p1.json`, 2026-10-19 08:15): the 12 wall corners read median 14.2 µm (4.4–31.5; the largest beside the arrows and label) against 11.0 drawn; the three fluid corners median 6.3 µm against 7.0.

**Possible solutions.** (1) A junction-corner radius drawn and modelled apart from the head fillet (rc7.x); (2) raise the Dinh preset's fillet to 8.4 µm (the fit's limit; it would also round the fluid corners to 8.4); (3) disclose (rc7).

**Status.** rc7: drawn (option `meander.junction_radius_um`, `Spec.r_junction`; `layout.assemble` rounds the two wall corners of every unit's stem where it leaves the trunk, mirrored), 11 µm in both Dinh presets (M); drawing only, since c_ent is an external input; LIMITATIONS 27. Every corner of the drawn head is rounded.

## 2. Device 2: Gifford et al. 2014, Fig. 4 (rc6 package: `gifford2014_fig4_design`, then the 51 µm start; rc7 package: `gifford2014_fig4_aligned`, reading B, with `gifford2014_fig4_design` as the reference and the rejected 51 µm start, `gifford2014_fig4_start51`, as a figure; RC8 package: `gifford2014_fig4_aligned`, the nested reading, with the reference and the rejected readings B (`gifford2014_fig4_pathorder`) and C, and the folded-head figure)

The rc6 package reproduced rc6 exactly:

- the DXF is byte-identical (`09b471f9…`), exit 0;
- every value of `gifford_results_rc6.json` is reproduced;
- the figures' cyan/amber split is the trace's own (88 / 106 edges).

The analyses below use the rc6 registration unchanged. The tools (`tools/gifford4_leg2.py`, `tools/gifford4_post_along.py`) read the paper image, the package's photomask and the rc6 record. They are maintenance tools, not part of the package.

### GIF4-1. Leg 2's side channels are about 18–19 µm wider per side on the image than the mask draws (decided 17:28: aligned; RC8: the nested reading, leg 2 open)

**What the image shows.**

- Across leg 2 (amber, not claimed), dark-line centres on 50 profiles read side channels of 124.5 and 125.0 µm (walls ± 233 µm from leg 2's axis).
- The mask draws 104.9–106.2 µm across the frame (walls ± 213–214 µm; the record rounds this to 104.5).
- The difference is a width, not a placement: both walls lie outward by the same ~20 µm, while leg 2's pills and centre line sit on the image's. A registration or spacing error would move one wall in and the other out.
- The trace record's own nearest-outline offsets of leg 2's walls agree: +17.5 to +18.2 µm per wall (`measured.leg2_side_walls`).
- Read the same way, leg 4 matches the print: 300.4 / 299.7 µm sides (printed final width 300) and a 100.7 µm centre (printed 100). So the method is not biased by 20 µm.

**What the printed design allows.**

- **Gifford's schedule.** His design is a constant 100 µm centre, a constant f_gap, and his duct law, eqn (5). Its side channels follow r(w_s(i)) = 2 r(w_c) q^i / (1 − q^i), with q = 1 − 2 f_gap (two gaps per row). That is Gifford's eqn (4) at equal pressure drop. opencif's exact closure gives the same widths to 0.1 µm, tested at seven rows.
- **The fraction does not depend on f_gap or N.** w_s depends on i only through q^i, so the fraction i/N at which a width is reached is ln Q(w)/ln Q(300 µm), with Q(w) = r(w)/(2 r(w_c) + r(w)).
  - 124.7 µm is reached at 57.5 % of the rows (57.3–57.9 % for depths 77–83 µm, the caption's 80 ± 3).
  - The mask's leg 2 at the inset sits at 48.6 %. Its widest row anywhere, its last, is at 50.4 % and is 108.4 µm wide.
  - So the image's leg 2 is wider than any row the drawn leg 2 holds.
- **No layout fixes it.** Leg 4 is the last leg: its array ends in the inset, with 300 µm sides and three leads to the outlets. Leg 2 meets the inset at leg 4's end. Then, with legs 1–3 of L rows plus a extra rows on leg 1, and leg 4 of L − s rows, a 124.7 µm side channel on leg 2 at the inset needs a − s ≈ 0.71 L. That is a first leg about 1.7 times the others, about 70 mm on a 47 mm chip.
- **The schematic makes it worse.** Fig. 4's schematic (not drawn to scale) puts leg 4's array end about a quarter of the way down the legs. That moves leg 2's rows at the inset earlier, so narrower, not wider.

**Rejected causes.**

- A modest start width. rc6's as-fabricated photomasks, measured at the inset, draw leg 2's sides 105.6 µm (8 µm start) and 106.6–107.9 µm (20 µm start). (Corrected 15:40 on 2026-10-18: an earlier line here gave 102.8 / 104.4 µm, from the planner run without the start side flow S0 that the pipeline sets.)
- The design depth. 124.7 µm would sit at about half the rows only for a design depth near 200 µm (48.7 %); Gifford's devices are 80–125 µm deep.
- f_gap and N (they drop out).
- The gap G (G does not enter the schedule; GIF4-3).
- The measuring method (leg 4 is right).

**The one mechanism inside Gifford's published procedure.** Gifford's eqn (4) runs from a chosen first width: "the desired value for ws(i = 1), which may be made close to or equal to zero, depending on the device production capabilities available" (p. 4498). Every later width inherits it. opencif implements this rule as `start: {style: floor, side_width_um: …}`; the design preset starts sealed (his Fig. 2 convention, w_s(0) = 0).

**Test variant, built 15:31 on 2026-10-18.** This is analysis only, not a package: `persist/devpkg/g4variant/`, the design preset plus `start: {style: floor, side_width_um: 52}`, rc6 code.

- **Where it matches.**
  - Leg 2's walls land on the image's: offsets −0.2 to −2.1 µm, against +17.5 to +18.2 µm for the design. So about 50 µm would centre them.
  - Leg 4 still matches the print.
- **What it costs.**
  - The design assumes a side flow S0 = 0.254 of the inlet flow per side at gap 1, so 34 % of the inlet flow is already in the side channels. A device whose side channels start sealed must pull that through its first gaps.
  - The rc6 network check on the drawn test mask: the first 37 rows realise cutoffs up to 11.1 × D (57.5 µm at row 1, so every bead passes there), and those rows carry 39 % of the filtrate.
  - N 2931 rows against 3546; the chip is about 40.4 mm long against 47.4.
- **Assessment against the established physics (A. H. Routt's question, 15:48).** The network solve is honest, but the design is not a CIF design by the frozen rules.
  - **The design rests on a fictitious flow.** Appendix A (frozen, "The normalised march"): with a sealed inlet σ0 = 0, and "The design code's other start mode prescribes a printable first side width; its σ0 > 0 is a fictitious inlet side flow, reported as a mass defect 2σ0". Here 2σ0 = 0.508 of the inlet flow.
  - **The drawn device breaks the constant-cutoff target** on which App. B's design map rests: 37 rows at up to 11.1 × D, carrying 39 % of the filtrate.
  - **What rc6 permits.** It allows a floor start only as a declared departure (Gifford's rule for fabrication floors, where small floors such as 8–20 µm move the entrance little).
  - **Verdict.** A 52 µm start is a physically computed but non-CIF design, not merely a departure from Gifford's text. Rejected.

**The layout, re-examined (15:55).**

- **The paper's arrows do not show leg 2's direction.** The two arrows in the inset (dark blobs of about 760 px) sit at mask y = 1136 and 2138 µm. That is 500 and 502 µm from leg 4's axis (1636), each about 92 µm outside one of leg 4's walls: a symmetric pair that most likely both mark leg 4. The arrow between the legs is also equidistant from leg 2's wall (85 µm), so it could mark either leg. The arrows are therefore neutral on the neighbour's direction.
- **The schematic's topology is unambiguous.** Fig. 4's schematic (enlarged, `gif5-001.png` columns 0–200) draws the inner loop (inlet → leg 1 down → U-turn → leg 2 up) inside the outer loop (top U-turn → leg 3 down → U-turn → leg 4 up to the outlets). That is nested, so leg 4's neighbour is leg 2, as the caption says.
- **The schematic's inset box sits about a quarter of the way down the legs.** Leg 4's leads (three channels, no posts) run from the outlets down to the box. Under the nested reading that puts leg 2 at the inset about 25 % before its end: row ≈ 1564, w_s ≈ 95.5 µm, which widens the gap to the image's 124.7.
- **A sequential layout fits.** If leg 4's neighbour were leg 3 (legs side by side in path order: 1, 2, 3, 4, a four-leg serpentine with the ports at the inlet end), the printed design with the sealed start reproduces the image.
  - With equal legs of L rows (leg 1 nine rows short at the inlet, as rc6 draws it) and leg 4's array ending s rows before the top: 4L − s = 3555, and leg 3 at the inset is at row ≈ (L − 9) + L + s + 3 = 2040.
  - So L ≈ 934 rows (42.8 mm legs) and s ≈ 179 rows (8.2 mm of leads before leg 4's array begins, from the outlet end).
  - The predicted s/L = 0.19 agrees with the schematic's box at about a quarter of the way down (not to scale). Under the nested reading the box position makes the mismatch worse.
  - The cost: it contradicts the schematic's drawn nesting and the caption's "second and fourth" (it would be "third and fourth").
  - rc6 cannot draw it: `_plan_paperclip` requires nested leg axes (0, +, −, +, …) and the serpentine requires an odd leg count.

**Possible solutions** (revised 15:55 on 2026-10-18).

| option | physics (App. A/B, glossary, rc6) | printed values | micrograph | schematic, caption | code | verdict |
|---|---|---|---|---|---|---|
| A. The printed design, nested (the package) | strict | all kept | leg 4 on; leg 2 20 µm inside per side | kept | rc6 | keep as the package; disclose |
| B. [Superseded: 17:28 (full alignment), then 00:44 (INV-2): reading B with the 8 µm floor start is rc7's Fig. 4 device; C (now 51 µm) is rejected by the beads; Fig. 6 is start40 with G 17.1 µm] The printed design, sequential legs, leg 4's array ending about 8 mm before the top | strict | all kept | leg 4 on; neighbour on (s read from its width through the physics: consistency) | topology and leg number contradicted; box position agrees | rc7 layout option (no physics change) | the only strict reading that reproduces the image; **weakened by device 3** (Fig. 6 is nested and shows a start offset): dropped as the rc7 candidate unless wanted |
| C. [Superseded: 17:28 (full alignment), then 00:44 (INV-2): reading B with the 8 µm floor start is rc7's Fig. 4 device; C (now 51 µm) is rejected by the beads; Fig. 6 is start40 with G 17.1 µm] 52 µm start | not a CIF design: fictitious inlet flow 2σ0 = 0.51 (App. A); 37 rows up to 11 × D | start "close to zero" contradicted | on (fitted on this leg: no test in Fig. 4) | kept | rc6 | rejected as opencif's design; **reconsidered (17:15) as the account of what Gifford built**: the mechanism is tested on Fig. 6 (GIF6-1); offer as a labelled as-depicted reproduction |
| D. A different f over legs 1–2 | the closure holds, but there are two cutoffs | the single printed f_gap contradicted | on (fitted) | kept | new code | rejected |
| E. Another f, N, depth, gap or centre width | strict | changed | cannot fix: the row fraction does not depend on f or N, depth would need about 200 µm, the gap does not enter, the centre reads 100–102 µm | — | — | ruled out |
| F. The as-fabricated 8 / 20 µm starts | small declared entrance error | kept | leg 2 +0 / +1.7 µm only | kept | rc6 | ruled out |

**Option B tested (16:10–16:17 on 2026-10-18; analysis only, `persist/devpkg/g4optionB/`).**

- **The setup.** The printed design preset is unchanged except for two layout values: `reproduction.leg_Y: [0, 818, 1636, 2454]` (legs side by side in path order) and `pc_last_short: 8200`. The code is the rc6 code with one check relaxed in a scratch copy: `_plan_paperclip`'s nested-order test. Loops take their own turning direction, and nothing in the physics changed.
- **The build.** It passes exactly as the nested design does:
  - exit 0, REPRODUCED; hard checks PASS; validity envelope INSIDE;
  - DRC at the floor: nothing lost; one fluid component; 7090/7090 holes; no part closer than the floor; ports separated;
  - Q1 and Q2 hold. Every gap realises D within ±1 % (5.1996–5.2067 µm, row 1 at 1.000 D), 100 % of F within 1 %, no declared rows, F 0.906475.
  - The junction band covers 10 rows at the last loop, up to 1.7 % of D (the nested design: 4 rows, up to 1.8 %).
  - Legs of 924 / 933 / 933 / 757 rows; chip 47.8 × 5.0 mm (nested: 47.4 × 5.0).
- **The trace.** This is the rc6 trace code, unchanged. A record copy gives the neighbour's axis (leg 3) in the slot the code reads as "leg 2". The registration is identical to the package's: V1 −0.09 / 0.04 px, median 1.21 µm, nulls 4.11 / 4.82 µm, leg spacing 817.9 µm on the image against 818 drawn.
  - The neighbour's walls: −0.78, −0.63, +1.09 and −1.34 µm (mean signed nearest-outline offsets by wall), against +17.5 to +18.2 µm for the nested design.
  - Leg 4 is unchanged: centre 100.4 µm; sides 299.9 / 298.7 µm.
- **What it shows, and what it doesn't.**
  - Under strict physics and every printed value, a sequential layout exists that passes every rc6 check and puts the micrograph's neighbouring leg on the image's walls. (This was the sealed form, at the design-level 4.5 µm floor, whose 4.86 µm first side channel is 16.5 : 1 deep against the 10 : 1 limitation Gifford prints, p. 4503; rc7's reading B takes his 8 µm floor instead, with rows 1–3 above D: INV-2.)
  - The leads' length (s = 8.2 mm) was read from that leg's width, so the wall match is consistency, not a test.
  - Its independent support is weak: the schematic's inset box about a quarter of the way down the legs, against a predicted 0.19.
  - Its cost: the schematic's drawn nesting and the caption's "second and fourth" would both be wrong.

**Recommendation (16:20, before device 3).** Keep A as the package for now. Carry B into rc7 as a layout option: legs in path order, ports at the inlet end, no physics change.
- If B is adopted, the Fig. 4 package would carry both A and B, traced. Which is primary is decided with the evidence of device 3. Decide which is primary only after device 3. Gifford's Fig. 6 micrograph shows four legs whose order rc6 checked against the nested layout, and Fig. 6 has its own printed inconsistency (final side width 140 µm in the caption, 144 µm in the body). That is independent evidence on how Gifford's lab laid out and sized these devices.

**Status (16:20).** open: A shipped. B tested, and it passes every check and matches the image. B is the rc7 candidate, and which of A and B is primary waits on device 3. C and D are rejected.

[Superseded: 17:28 (full alignment), then 00:44 (INV-2): reading B with the 8 µm floor start is rc7's Fig. 4 device; C (now 51 µm) is rejected by the beads; Fig. 6 is start40 with G 17.1 µm] **Revised with device 3 (17:15 on 2026-10-18).** Fig. 6 is the independent evidence the recommendation waited for (GIF6-1).

- **Fig. 6 is nested.** Its four visible legs all flow the same way (four arrows, one per leg), as the nested layout requires; in a serpentine, neighbouring legs alternate. Where the paper shows the direction, the layout is the schematic's.
- **Fig. 6 shows a start offset, tested in part (GIF6-2).** One start width (40 µm), fitted on one leg with the analytic schedule, predicts the other legs' widths. Drawn by rc6 and traced, every leg end lies within 1.5 µm of the image (the sealed design: 5.4–10.0 µm narrow), and every leg's registration improves. The strict alternative there (23 mm of extra rows) has no support.
- **What that does to Fig. 4's options.**
  - *C (nested, about 52 µm start)* is the mechanism observed and tested in part (GIF6-2) in the same paper's other device, and it keeps the schematic's nesting and the caption's "second and fourth". On Fig. 4 itself it remains a one-leg fit (no test). Its cost is unchanged: a declared departure (App. A's mass defect 2σ0 = 0.51), 37 rows up to 11.1 × D carrying 39 % of the filtrate.
  - *B (sequential)* stays strict, but it now needs the schematic, the caption and the pattern of the paper's other device all to be wrong about this one.
- **The physics verdict is unchanged.** C is not a CIF design by the frozen rules; it is Gifford's departure, which rc6 models exactly as a declared `start: floor`. The earlier "rejected" was a verdict on C as opencif's design. As an account of what Gifford built, the evidence now favours it: the mismatch is Gifford's (his start rule, unreported), not opencif's.

**Recommendation (17:15).**

1. Keep A as the package's device for Figs. 4 and 6: strict physics and the printed values (A. H. Routt's rule of 15:48).
2. Attribute the leg-width discrepancies of both micrographs to Gifford's start rule (C for Fig. 4, start40 for Fig. 6), disclosed in the packages and the dissertation.
3. Drop B as the rc7 layout candidate, unless A. H. Routt wants it kept as a strict alternative.
4. Optional, for rc7: as-depicted presets `gifford2014_fig4_start52` and `gifford2014_fig6_start40` (declared departures, their cost reported), traced as secondary, labelled reproductions in the two packages, so the micrographs can be shown matched. rc6's as-fabricated presets (assumed 8 / 20 µm and 15 µm starts) could take the measured starts instead.

**Status.** [Superseded: 17:28 (full alignment), then 00:44 (INV-2): reading B with the 8 µm floor start is rc7's Fig. 4 device; C (now 51 µm) is rejected by the beads; Fig. 6 is start40 with G 17.1 µm] open: a decision for A. H. Routt at rc7 (items 3 and 4).

**Decided (17:28): full alignment.** The reproduction of Fig. 4 is C with the measured gap: the printed design (f_gap 3.34 × 10⁻⁴, w_c 100, depth 80, final width 300 µm, nested) with a 52 µm start (**M**, fitted on leg 2) and G 17.2 / g_p 28.6 µm (**M**, GIF4-3). B is dropped; A stays as the design-level reference. Tested 17:24–17:29 (analysis only, `persist/devpkg/g4complete/`): exit 0, REPRODUCED with declared departures (floor start); leg 2's walls −0.3 to −2.1 µm from the image's (A: +17.5 to +18.2); leg 4 unchanged (centre 100.4, sides 299.9 / 298.7 µm); the posts' flanks along the flow +0.26 µm median (A: +1.67), the lead dividers' ends +0.1 to +0.5 µm (A: −1.3 to −1.4), the claim median 0.63 µm and p90 1.41 µm at the ± 3.5 px window (A: 1.09 / 2.51); the rc6 ± 8 px window still reads p90 16.6 µm (GIF4-2). Its entrance, reported: rows 1–45 more than 1 % above D, row 1 at 10.0 × D, 39 % of the filtrate. rc7: preset `gifford2014_fig4_aligned`, traced; the package's primary device. **Superseded (2026-10-19 00:44, INV-2):** C is rejected by the inset's beads; `gifford2014_fig4_aligned` is reading B, and C is kept as `provenance/rc7_traces/variants/gifford2014_fig4_start51.yaml` (51 µm), shown as a figure.

**RC8 (2026-10-21): the nested reading; reading B rejected by the inset's taper direction.**

- **The test (V13).** In a CIF leg the side channels widen along the flow (Gifford's Fig. 1 caption: "The width of the
  side channels, ws(i), gradually increases with increasing gap row number, i"), so the sign of a leg's outer
  wall-to-wall width slope is its flow direction.  The neighbour's outer walls were located on the inset through each
  reading's own registration, at every 5 µm along the leg inside the micrograph (`provenance/rc7_traces/code/
  taper_sign.py`), and the image's width and the mask's were fitted against the position: nested (the neighbour is
  leg 2) PASS: the image's slope −0.00240 ± 0.00059 µm/µm (4.1 SE) against the mask's −0.00257, 130 of 130 eligible positions at a ±30 µm window; reading B (leg 3) FAIL: −0.00254 ± 0.00062 (4.1 SE) against +0.00301, 117 of 117; leg 4's own array part (the control) PASS: −0.00469 ± 0.00040 (11.6 SE) against −0.00800, 80 of 80 at ±12 µm on the nested mask and PASS on
  every other.  The image's neighbour narrows the same way in every reading (each map registers leg 4 on the same
  image); only the mask's sign changes, and B's is opposite.  First measured by Review 1 of rc7 (its validator V13,
  −0.0028 ± 0.0003 µm/µm), re-measured with opencif's own implementation and adopted as a sign test: the magnitude
  moves by about 15 % with the sample grid and is not judged; the standard error is the larger of a Newey-West one (Bartlett, lag 100 µm,
  beyond two row pitches) and the independent-sample one (neighbouring positions share pixels: the residuals' lag-1 correlation is 0.64–0.77, so the independent-sample SE understates the
  uncertainty: on the neighbour the Newey-West SE at lag 100 µm is 1.9–2.0 × the independent-sample one; on leg 4's
  controls it is 0.8 ×, and the recorded SE there is the independent-sample one; the RC8 code review, M1, N2); the wall window must exceed the largest wall offset (±30 µm on the neighbour:
  the nested reading's walls lie 16.3–19.6 µm outside the mask's; at ±12 µm none is found and the test is
  inconclusive: review 3, E6).  A row-spacing control (the image's row pitch against the mask's along the
  neighbour) was inconclusive (−0.0031 ± 0.0023), so the sign of the width slope is the evidence, not the pitch.
- **The decision (D1, D6).** The reproduction is the nested reading with the 8 µm start and the inset's measured gap,
  under the name `gifford2014_fig4_aligned`: every measured value reading B had (G 17.2 / g_p 28.6 µm, the legs
  818 µm apart, the registration) and the layout the evidence supports: the caption ("the second and fourth legs");
  the schematic's topology; the inset's two arrows, which point the same way (the weakest line: their placement,
  symmetric about leg 4, also allows both to mark leg 4); and the taper direction.  Leg 4's array ends 2 mm short
  (rc6's reading of the inset's leads).  Built at the RC8 code: exit 0, REPRODUCED with declared departures (floor
  start), N 3543, F 0.90650, rows 1–3 above D (row 1 at 1.97 D, 0.47 % of F), 99.5 % of F within 1 % of D.  On the
  inset: the claim (leg 4) median 0.61 µm, p90 1.40 µm, coverage 0.999, censored null ratios 8.66 / 6.46, chamfer null
  ratios 2.89 / 2.38 (reading B: 0.62 / 1.40).
- **What it does not reproduce (open).** Leg 2's side channels beside the inset read 17.5–18.1 µm per side wider on the
  image than the printed schedule draws them (the trace record's `measured.leg2_side_walls`, every wall class with at
  least 50 samples).  On the nested mask the rows that pass the inset are fixed: leg 2's rows 1730–1758 lie inside the
  micrograph, drawn 105.1–106.7 µm wide; the image's 122.6–124.8 µm (the drawn widths plus the walls' measured offsets; the
  width profiles read 124.5 / 125.0 µm) are what the schedule draws at rows 2008–2040, in
  leg 3, about 280 rows (12.7–13.0 mm of leg) further on; as a scale, 1.17 × the schedule's width at row 1744
  (`verification_here/review3_checks/leg2_check.py` in the RC8 reviews zip; review 3 found rows 1717–1764 against 1996–2045 with a wider window, the same in
  substance).  A registration or wall-finding error cannot be that large: the pitch and the legs' spacing read to
  0.1 %, leg 4's walls within one image pixel, and both of leg 2's walls lie outward by the same amount while its
  posts and centre line sit on the image's.  No reading with the printed single f and a start inside the beads'
  bracket explains it.  The dissertation states it beside the claim.
- **Reading B** (rc7's `gifford2014_fig4_aligned`) is the variant `provenance/rc7_traces/variants/gifford2014_fig4_pathorder.yaml`:
  built, traced and shown as a rejected reading.  Its mean wall offsets (−1.4 to +1.0 µm) agree only because leg 4's
  early end was fitted to the neighbour's width, and a mean over the window cannot see a slope.
- **Status.** RC8: decided (D1, D6); leg 2's width open, disclosed (LIMITATIONS 28, PRESET_PROVENANCE, DISSERTATION_NOTES
  § 3 item 7 and § 4, the preset's header, the Fig. 4 package's README).

### GIF4-2. The claim's p90 (18.2 µm) and its −12.2 µm class are the metric's search window, not geometry (rc7, trace record)

**Cause.** The nearest-outline window is ± 8 px = ± 22.3 µm at 2.79 µm/px. That is wider than the 20.5 µm gap between leg 4's posts. Along the flow, a sample whose own outline is faint finds the neighbouring post's outline.

**Evidence (`tools/gifford4_post_along.py`).** The tool reproduces the record exactly at ± 8 px (n 7054, median 1.21, p90 18.17 µm).

- On the posts' flanks (n 1700, the record's 356 + 1344), the offsets fall into two groups:
  - 1089 samples (64 %) within −3 to +6 µm;
  - 604 (36 %) at +15 to +24 µm, the next post.
- With ± 6 px (16.7 µm) or ± 3.5 px (9.8 µm):
  - the flanks read median +1.65 µm, p90 |d| 2.8 µm, none beyond 15 µm;
  - the lead dividers' ends (the record's channel_along_C, −12.2 µm, n 39) read −1.3 to −1.7 µm;
  - the whole claim reads median 1.04–1.09 µm, p90 2.49–2.51 µm, with coverage 48–54 %.

**What it settles.** The three rc6 texts on this tail are all superseded (adjudication § 5.4, edits A4b and A5): "≈ 4 µm long", "about 7.5 µm … on every post alike", and the README's mechanism.

**rc7.**

- A rule in the trace code: a panel's window must be shorter than the smallest gap between outlines along the sample's normal (G here), and the record carries the window, as it already does (`metric_search_px`).
- Retrace Gifford Fig. 4 at ± 6 px (or ± 3.5 px), which makes a new reference record.
- Check every other panel against the rule:
  - Dinh: ± 8.5 µm against G 19 µm, fine.
  - Gifford Fig. 6: ± 11.3 µm against 17.1 µm, fine.
  - **Xia: ± 14.4 µm against G 16 µm, too close.** Its p90s of 12–13 µm need the same test in its package.

The trace code is provenance, not hashed, so the code stamp does not change. The records and the docs that quote them do.

**Done (rc7).** The rule as implemented is stricter than "shorter than the smallest gap": a half-window at most half the gap (`panel_core.window_px`): Dinh 8.53 µm, Fig. 4 8.37 / 9.77 µm (aligned / design), Fig. 6 6.79 / 9.05 µm, Xia 7.19 µm (f 7.98 µm). Under it Fig. 6's 11.3 µm against 17.1 µm fails ("fine" above is the first rule's verdict), so rc7 uses 6.79 µm (`provenance/rc7_traces/README.md`).

### GIF4-3. The gap reads 17.2 µm along the flow on the inset; the preset uses 20.5 µm (decided 17:28: 17.2)

**The reading.** With the narrow window (GIF4-2), leg 4's posts read +1.65 µm outside the drawn posts at both flanks. That makes the posts 28.6 µm along the flow and the gaps 17.2 µm at the measured row pitch of 45.82 µm. The lead dividers' ends read solid-ward by 1.3–1.7 µm too. Gifford's Fig. 6 micrograph reads 17.1 µm (PRESET_PROVENANCE).

**Against it.**

- Across the flow the same image reads only +0.45 µm (posts) and ± 0.5 µm (walls).
- Gifford prints the gap as "~20 µm" (p. 4496).
- rc6 keeps 20.5 µm ("the record's measurement on Fig. 4"), treating the along-flow outlines as relief-biased.
- On a relief image a symmetric outline bias and a real difference cannot be separated.

**What the choice changes** (computed with the rc6 planner). At constant f, G does not enter the side-channel schedule. G 17.2 against 20.5 (g_p 28.6 against 25.3, the pitch fixed) changes no width by more than 0.001 µm, and N and F not at all. It changes:

- R_gap (+57 %);
- the drawn posts and gaps, and so the photomask;
- the gap diagnostics (Reynolds number, aspect ratio).

**Device 3's evidence (GIF6-3).** Fig. 6 reads the same gap, 17.1 µm, with the same anisotropy (+1.3 to +1.6 µm along the flow, +0.2 to +0.3 across), at 1.6 times the pixel size and nearly twice the depth. An outline bias fixed in pixels, or a slot shadow that grows with depth, would not give the same µm in both. This favours a real gap of about 17.2 µm in both devices, though it does not prove it.

**Device 4 (XIA-6).** Xia's fluorescence panels read the posts about 2 µm outside the drawn ones along the flow and 2–4 µm across, an isotropic outline bias, so Xia's printed gap stays. That neither supports nor contradicts this reading: Gifford's bright-field images show the excess along the flow only (+1.3 to +1.7 µm, against +0.2 to +0.45 µm across).

**Possible solutions.**

1. Keep 20.5 µm and disclose the reading (rc6's position).
2. Set G 17.2 µm, g_p 28.6 µm (**M**, both Gifford micrographs agree). This changes the photomask's posts, not its schedule.

**Status.** decided (17:28): set G 17.2 µm, g_p 28.6 µm (**M**) in the aligned preset (GIF4-1); the design-level preset keeps 20.5 µm as rc6's reference.

### GIF4-4. Layout values read on the same inset (note)

- **Leg spacing.** 818 µm, read on this inset; the generator's own rule gave 969 µm.
- **Leg 4's early end.** 2 mm (`pc_last_short`).

Both are read on the inset the trace uses, so where leg 2 lies there is consistency, not a test (LIMITATIONS 24 b). The port radii are assumed (the schematic is not to scale). No change.

### GIF4-5. Design-level start (note)

The sealed start draws a 4.86 µm first side channel against D = 5.20 µm, and the floor (4.5 µm) is chosen so that it passes. The build reports the aspect ratio of 16.5 (over 10). LIMITATIONS 26. This is not visible at the inset. No change: the as-fabricated presets exist for a fabricable start.

---

## 3. Device 3: Gifford et al. 2014, Fig. 6 (rc6 package: `gifford2014_fig6_design`; rc7 package: `gifford2014_fig6_aligned`, the 40 µm start, with `gifford2014_fig6_design` as the reference)

The rc6 package reproduced rc6 exactly:

- the DXF is byte-identical (`e47f5ca6…`), exit 3 (REPRODUCED outside the validity envelope, as rc6's record);
- every value of both rc6 records is reproduced: the four per-leg registrations (`gifford6_results_rc6.json`) and the one map with the turns (`gifford6_results_rc6_onemap.json`);
- the figures' cyan/amber split is the trace's own (one map 416 / 20 edges; each leg 106 / 6).

The analyses below use the rc6 records unchanged. `tools/gifford_row_fraction.py` uses rc6's duct law (`cifsched.rr`) and the record's widths; the test variant is analysis only (`persist/devpkg/g6variant/`: preset, `run_variant.py`, its output in the session scratchpad).

### GIF6-1. The side channels at the leg ends are wider on the image than drawn: by 9.6 µm at leg 1, falling to 5.5 µm at leg 7 (decided 17:28: aligned)

- *What the image shows.* rc6's record reads the side channels at each leg's end (dark-line centres; `measured.side_channels_at_leg_ends`): leg 1 58.5 / 59.3 µm, leg 3 87.1 / 87.5, leg 5 114.1 / 114.5, leg 7 145.5 (one side in frame). The two sides of a leg agree within 0.8 µm. The mask draws 49.3, 80.7, 108.7 and 139.9 µm. rc6 recorded this (LIMITATIONS 24 e) and named three possible causes: the assumed gap, an outline bias, or the device.
- *What the printed design allows.* Gifford's schedule (constant 125 µm centre, constant f_gap, sealed start) reaches a width w at the row fraction ln Q(w) / ln Q(w_end), with Q(w) = r(w) / (2 r(w_c) + r(w)); f_gap, N and G drop out. At the legs' ends (rows 861, 2601 and 4341 of 6080) it gives 49.4, 80.8 and 108.8 µm: the drawn widths to 0.1 µm. So the drawing is the printed design, and two of rc6's three causes are excluded:
  - **the gap** does not enter the schedule at constant f;
  - **a constant outline bias** cannot make an excess that depends on the leg: to put each leg on the schedule it would have to be 10.3 µm at leg 1, 7.1 at leg 3 and 5.6 at leg 5 (with the final width shifted by the same bias). Depths of 145–155 µm (the caption's ± 5) do not change this: the bias would still have to be 9.9–10.8 µm at leg 1 and 4.8–6.5 µm at leg 5.
- *What fits.* The excess falls along the device the way a **start offset** does: Gifford's own start rule, "the desired value for ws(i = 1), which may be made close to or equal to zero, depending on the device production capabilities available" (p. 4498), with every later width inheriting it (his recursion then follows the sealed schedule from the row where it reaches that width). One start width, fitted on leg 1 alone with the text's final width of 144 µm, is 40.0 µm. It predicts legs 3 and 5 at 86.7 and 113.6 µm against the image's 87.3 and 114.3 (−0.6 and −0.7 µm, a seventh of an image pixel), and leg 7 at 144 against 145.5. With the image's own final width (145.5 µm) the start is 39.7 µm and all three legs agree within 0.2 µm. With the caption's 140 µm, legs 3 and 5 miss by −1.6 and −3.0 µm. (`tools/gifford_row_fraction.py`, outside this package.)
- *Drawn by rc6 (analysis only, not this package).* A test build with rc6's existing start option (`start: {style: floor, side_width_um: 40}`) and F from the text's 144 µm, traced by the same rc6 code (`persist/devpkg/g6variant/`):
  - the drawn side channels at the leg ends are 58.8, 86.7, 113.6 and 143.9 µm against the image's 58.5 / 59.3, 87.1 / 87.5, 114.1 / 114.5 and 145.5 µm: −1.5 to +0.3 µm (the design: −5.4 to −10.0 µm);
  - every leg's registration improves: medians 0.87–1.35 µm (design 1.35–1.88), p90s 2.2–2.5 µm (4.3–8.8), the walls' mean signed offsets +0.1 to +1.7 µm (+4.5 to +9.0); one map 1.11 / 2.86 µm (1.47 / 5.93); the turns stay on (median 1.6 µm, the design 1.4);
  - the start was fitted on leg 1 with the analytic schedule; legs 3 and 5 test it in part (the text's 144 µm, which sets leg 7's end through F, was preferred over the caption's 140 µm with legs 3 and 5 in view: GIF6-2), and they pass; leg 7 is not a test.
- *What it means physically.* A 40 µm start is not a sealed CIF start. App. A (frozen): "its σ0 > 0 is a fictitious inlet side flow, reported as a mass defect 2σ0"; here 2σ0 = 0.11 of the inlet flow. The drawn device, sealed at its entrance, has to pull that flow through its first gaps, which then extract more than the design's f: in the test build rows 1–32 realise cutoffs more than 1 % above D, up to 11.0 × D at row 1 (39.9 µm, larger than any platelet), and those rows carry 14 % of the filtrate; the device needs 5,783 rows instead of 6,080 (Fig. 4's 52 µm start: 37 rows, up to 11.1 × D, 39 % of the filtrate). It is a declared departure that rc6 already models (`start: floor`), not a change to the physics.
- *Why it is likely what Gifford built.* The sealed design's first side channel is 4.69 µm wide in a 150 µm deep device (aspect ratio 32): it cannot be made, and Gifford's rule chooses the first width by "the device production capabilities available". rc6's as-fabricated preset assumed 15 µm without a measurement; this micrograph measures about 40 µm (aspect ratio 3.75).
- *The strict alternative.* The same widths follow, under strict physics, if about 456 real rows (23 mm of array) lay upstream of leg 1's end in addition to the seven equal legs. Leg 1 can be at most about 1.4 mm longer than the others inside the turn around it, and the schematic shows no extra fold, so this reading has no support.
- *Conclusion.* The micrograph's device was most likely built with Gifford's start rule at about 40 µm: a departure from strict CIF at its entrance, which his paper allows but does not report. The package keeps the printed sealed design, as rc6 does; the discrepancy is noted, not changed.

**Possible solutions.**

| option | physics (App. A/B, glossary, rc6) | printed values | micrograph | schematic, caption | code | verdict |
|---|---|---|---|---|---|---|
| A6. The printed design, sealed start (the package) | strict | all kept (the caption's 140 µm) | leg ends 5.4–10.0 µm narrow | kept | rc6 | keep as the package's device; disclose |
| start40. Gifford's start rule: a 40 µm start, final width 144 µm (the text) | a declared departure: App. A's fictitious inlet flow, 2σ0 = 0.11; rows 1–32 up to 11.0 × D, 14 % of the filtrate | f_gap, w_c and depth kept; the final width the text's; the start **D, fitted** (one value, on leg 1's widths, M) | every leg end within 1.5 µm (legs 3 and 5 a partial test; leg 7 set by F); every leg's median and p90 improve | kept | rc6 (`start: floor`) | the account of the device that the image supports; Gifford's departure, not opencif's: note it, and offer it as a labelled as-depicted reproduction |
| H6. Strict physics, with about 456 more real rows upstream of leg 1's end | strict | kept | on (the same widths) | contradicted: the schematic's seven legs show no extra fold, and leg 1 can be at most ~1.4 mm longer | new layout | no support |
| A constant outline bias | — | — | fails: it would have to be 10.3 / 7.1 / 5.6 µm at legs 1 / 3 / 5 | — | — | ruled out |
| Another f, N, G or depth | strict | changed | cannot: f, N and G drop out of the row fraction; depths 145–155 µm leave it | — | — | ruled out |

**An observation, not interpreted.** Both fitted starts put row 1's realised cutoff at 11.0 × D (Fig. 6: 11.01; Fig. 4: 11.05; test builds at G 20.5 µm: the final records give 9.90 D for Fig. 6 aligned and 9.89 D for the start51 variant), although the start widths (40 / 52 µm), depths (150 / 80 µm), mass defects (0.11 / 0.51) and cutoffs relative to w_c (0.32 / 0.58) all differ. With two devices it may be chance.

**Status.** [Superseded: 17:28 (full alignment), then 00:44 (INV-2): reading B with the 8 µm floor start is rc7's Fig. 4 device; C (now 51 µm) is rejected by the beads; Fig. 6 is start40 with G 17.1 µm] note: A6 shipped. start40 tested in part (analysis only, `persist/devpkg/g6variant/`). Whether rc7 adds as-depicted Gifford presets (start40 here, C in GIF4-1) as secondary, labelled, traced reproductions is a decision for A. H. Routt.

**Decided (17:28): full alignment.** The reproduction of Fig. 6 is start40 with the measured gap: the printed design with the text's final width 144 µm (GIF6-2), a 40 µm start (**D, fitted** to leg 1's widths, M; the package's measure step reads 40.1 µm from leg 1 and predicts legs 3 and 5, a partial test) and G 17.1 / g_p 34.3 µm (**M**, GIF6-3). A6 stays as the design-level reference. Tested 17:24–17:34 (analysis only, `persist/devpkg/g6complete/`): exit 3 (the same envelope verdict as A6), REPRODUCED with declared departures (floor start); every leg end within 1.5 µm of the image (A6: 5.4–10.0 µm narrow); each leg's median 0.51–0.64 µm, p90 1.55–1.74 µm (A6: 1.35–1.88 / 4.3–8.8); the posts' flanks along the flow −0.22 to +0.15 µm (A6: +1.2 to +1.6); one map 0.83 / 2.25 µm (final rc7 trace 0.82 / 2.19; A6: 1.47 / 5.93); the turns' walls a median 1.6 µm from the mask's; leg 7's drawn end 2 rows beyond the image's (rc6 places it in steps of about 6 rows: GIF6-4). Its entrance, reported: rows 1–40 more than 1 % above D, row 1 at 9.9 × D, 14.5 % of the filtrate. The package (version 2) builds it as its device.

### GIF6-2. The final side width: 140 µm in the caption, 144 µm in the text (decided 17:28: 144)

The preset uses the caption's (F 0.71768; the text's gives 0.72923). The image reads 145.5 µm at leg 7's end (one side), a third of a pixel from 144. Under GIF6-1's start offset the image's widths fit 144 (residuals ≤ 0.7 µm) and fit 140 worse (up to 3.0 µm). On the other side, the caption's "~2.5×" flow ratio is F / (1 − F) = 2.54 at 140 and 2.69 at 144. The platelet counts of the typical run (384 → 1039 central, 93 side; yield 83.2 %) balance at F = 0.69 (2.25×), below both: the delivered split also depends on the outlets and tubing, which the paper does not give. **Decided (17:28): 144 µm (the text, supported by the micrograph); the caption's 140 µm is stated in the dissertation as one of the paper's caption-text conflicts.**

### GIF6-3. The gap reads 17.1 µm along the flow; the preset assumes 20.5 µm (Fig. 4's) (decided 17:28: 17.1)

rc6's measurement reads 17.08 µm (16.9–17.2), and the end checks 16.7–17.9 µm. In the trace record the posts' flanks along the flow lie +1.3 to +1.6 µm outside the drawn posts, the posts' tips across the flow only +0.2 to +0.3 µm. Fig. 4's inset reads the same: 17.2 µm, +1.65 along against +0.45 across (GIF4-3). Fig. 6 has 1.6 times the pixel size (4.52 against 2.79 µm) and nearly twice the depth (150 against 80 µm). An outline bias fixed in pixels would differ by 1.6 times in µm, and a narrow slot's shadowing would grow with its depth; neither is seen. That favours a real gap of about 17.2 µm in both devices, or a cause that depends on neither. At constant f the gap does not enter the schedule. Added to GIF4-3 as evidence. **Decided (17:28): G 17.1 µm, g_p 34.3 µm (M) in the aligned preset.**

### GIF6-4. Leg 7's drawn array ends beyond the image's: 4 rows (205 µm) on the design preset, 2 rows (103 µm) on the aligned one (rc7, layout)

The image's leg 7 ends 5 rows (257 µm) before the other legs. The generator ends every middle leg at a whole row and gives the last leg the remainder, so `pc_last_short` moves leg 7's end only in steps of about 6 rows: 257 µm gives 1 row (LIMITATIONS 25). rc6's one-map claim stops 6 rows before leg 7's drawn end for this reason (amber in the figures). *rc7 solution (proposed):* let the first leg take the remainder (its start at the inlet moves by whole rows), so the last leg's end can be placed exactly; a layout change only, no row's physics changes. **Status: not done in rc7** (the layout still gives the last leg the remainder, LIMITATIONS 25); rc7's trace code instead derives the one-map claim's end from the drawn mask (`LEG7_IMAGE_SHORT_ROWS`, `LEG7_MARGIN_ROWS`; the record's `measured.leg7_claim_end`; package fix 1).

### GIF6-5. Legs 2, 4 and 6 are placed by assuming semicircular turns (note)

Only 8–23° of each turn is in the frame. The turns' walls then lie a median 1.4 µm from the mask's (consistency: the axes were read on those walls). The independent result is along the flow, which nothing read on the turns sets: the mask's turn centre lies 5.7 µm downstream of the image's. No change.

### GIF6-6. The design-level start (note)

The sealed start draws a 4.69 µm first side channel at 150 µm depth (aspect ratio 32), with the floor chosen so that it passes (LIMITATIONS 26). Not visible on the micrograph, but GIF6-1 measures what the device used instead. No change to the design preset; the as-fabricated preset's assumed 15 µm start could become the measured ~40 µm at rc7 (not done in rc7: `gifford2014_fig6_aligned` carries the 40 µm start instead).

### GIF6-7. Length (note)

The text gives "~40 cm, in total". The design draws 35.9 cm of array path (rows and turns; `drawn_parameters.array_length_mm` 359.2) on a 53.8 × 8.3 mm footprint (4.5 cm², within the text's "below 6 cm2"). The 40 µm start of GIF6-1 needs 5,783 rows instead of 6,080: 34.4 cm. Both are "~40 cm" only loosely; not a test. No change.

## 4. Device 4: Xia et al. 2016, RS-7 (rc6 package: `xia2016_rs7_depicted_fig3` and `xia2016_rs7_depicted`; rc7 package: `xia2016_rs7_aligned`, the parallel-plate taper, with `xia2016_rs7_depicted` as the reference)

The rc6 package reproduced rc6 exactly, for both photomasks:

- the aligned device, rc6's variant `xia2016_rs7_depicted_fig3` (D 6.5705 µm and F 0.94072 read off Fig. 3a and panel g): DXF byte-identical (`280a04ed…`), exit 3;
- the printed reference, rc6's preset `xia2016_rs7_depicted` (D 6.8644 µm from the printed f_ref, F 0.938 from the printed ~15:1): DXF byte-identical (`0161ab00…`), exit 3;
- every value of rc6's twelve trace records: the whole device on Fig. 3a (V11, V12), panel g, panels c, d, e, f, for each photomask.

**The decision this package carries (A. H. Routt, 17:25–17:33: full alignment).** The open choice "the shipped preset or the fig3 variant" is settled by the aim: the fig3 variant is the package's device (it aligns Fig. 3a's whole device: far loops +0.09 mm, outlet ring −0.19 mm, against +5.69 / +5.41 mm for the printed preset; and panel g's last side channel). Its D and F are fitted to the figures (**D, fitted**: F to panel g's last side channel, D to Fig. 3a's leg fill); the text's values (printed f_ref → D 6.8644 µm; "approximately 15:1" → F 0.938) are the reference. For the dissertation: the figures need a 4.3 % smaller cutoff and a 15.9:1 ratio, which the text's "approximately 15:1" allows and its f_ref, through the frozen cutoff law, does not quite. *(rc7 and RC8: with XIA-1's parallel-plate taper the same fit gives D 6.4159 µm, 6.5 % below the printed f_ref's 6.8644 µm, and F 0.94072 as before; `presets/xia2016_rs7_aligned.yaml`. The 4.3 % is rc6's fig3 variant's, D 6.5705 µm.)*

**What rc6 cannot align, and why.** Three of the figure's features are off on both photomasks: the centre channel's taper (XIA-1), panel f's side channels (XIA-3), and Fig. 1e(ii)'s loop-1 widths (XIA-5). XIA-1's cause is found and needs an rc7 taper option; XIA-3 is open.

### XIA-1. The centre channel narrows faster than rc6's constant-shear taper (rc7: a taper option)

**What the figure shows.** The centre width read on the panels (rc6 records, 1.798 µm/px): 299.9 µm at leg 1's first posts (panel c), 219.8 µm at leg 3's first posts (panel d), 186.0 µm half way down leg 3 (panel e, placed by Fig. 3a's box), 150.3 µm at leg 5's first posts (panel d), 151.5 µm on legs 5 and 7 (panels e, f), 148.3 µm at the outlet (panel g). The aligned photomask draws 300, 230.9, 205.5, 150, 150, 150, 150 µm. So leg 3 is drawn 11–20 µm too wide; the rest matches. Panel e's two registrations disagree across by 3.7 px (V1 fails) for this reason: its post columns are 19 µm closer together than drawn.

**What Xia prints.** The principle and the ends, not the profile: "determined numerically (using formulas for rectangular channels based on Yang et al.) in order to maintain a relatively constant amount of (maximum) shear experienced along each leg of the device" (p. 4); 300 µm at the inlet, "gradually decreases until it reaches 150 µm" (Fig. 1d(i)); Fig. 2b's maximum shear roughly constant over the first half of the device. rc6 draws the principle with the exact Stokes maximum wall shear of the centre channel, γ = C(α) Q / (a b²) (spec v5 S16-b), which is **A** in its form (LIMITATIONS 16: no user profile).

**The cause (analysis, 18:00–18:15; `persist/devpkg/xiataper/`).** A taper w_c ∝ Q_c^n, held with the other two readings of the figure (N 8106 rows, Fig. 3a's legs filled; F from panel g's 624.6 µm), solved with rc6's own planner in a scratch copy that only adds the exponent:

| taper | leg 3's first posts (image 219.8 µm) | half way down leg 3 (image 186.0 µm) | leg 5's first posts (image 150.3) | D (fills Fig. 3a) | F |
|---|---|---|---|---|---|
| rc6's exact constant shear (the aligned device) | 230.9 | 205.5 | 150 | 6.5705 µm | 0.94072 |
| n = 0.90 | 225.5 | 197.3 | 150 | 6.4948 | 0.94072 |
| n = 0.977 (fitted on leg 3's first posts) | 219.8 | 188.7 (test: +2.7) | 150 | 6.4325 | 0.94071 |
| **n = 1: w_c ∝ Q_c** | **218.1** | **186.1** | **150** | **6.4159** | **0.94071** |

n = 1 is not a fit: it is constant maximum shear by the parallel-plate formula, γ = 6 Q / (w_c h²), the coefficient C(α) held at its wide-channel limit. It puts both of leg 3's panels within one image pixel (−1.7 and +0.1 µm), where rc6's exact form is 11 and 20 µm off. **But it departs from spec v5 S16-b, and Xia's text does not state it** (foundation checks, 2026-10-19): Xia's "formulas for rectangular channels based on Yang et al." are rectangular-duct formulas, whose Eqn 3 the glossary identifies as the resistance series (G:2036–2038); the text does not state the shear coefficient. Spec v5 holds constant shear with C(α) (S:456) and reports "max gamma_dot via C(alpha) [NEVER 6]" (S:606). So the figures read a form spec v5 excludes, which Xia's text neither states nor rules out: a reproduce-mode departure, labelled as a reading of the figures, not as Xia's method. The record reports the true peak shear via C(α): 1.34–1.57 × the parallel-plate value, rising 17.4 % along the taper (`validation.gates.shear_peak`).

**Physics.** Any taper keeps every row at D (spec v5 S1, "exact at every aspect ratio … ∀ 𝒲"; App. A's march and closure hold for any w_c(i)), so the change is the taper's form, not the physics. Its cost: D falls to 6.42 µm to keep Fig. 3a's legs filled (more extraction per row on legs 1–4 at the same D means fewer rows; the fill asks for the same rows).

**Possible solutions.**

1. **(Done, rc7)** A taper option for this form: `centre_channel.taper: constant_shear` with `shear_formula: exact | parallel_plates` (default `exact`, so every rc6 photomask stays byte-identical; `parallel_plates` in reproduce mode only, a departure from spec v5 S16-b). The aligned Xia preset takes `parallel_plates` (**A**, chosen with panels d and e in view, which are then consistency, not a test), with D and F refitted (6.4159 µm, 0.94072), built and traced (`review7/xia1/`): panel e's registrations now agree (V1 across 0.94 px; rc6 −3.73), panel d's p90 6.3 µm (9.96; the final rc7 trace, at the window rule's 7.19 µm: 4.76); panel f unchanged (INV-4).
2. A user taper profile (w_c at measured positions, interpolated): more general, but it fits the four panels rather than testing them. Not preferred while one printed principle fits.
3. Keep rc6's form and disclose (the state of this package).

**Status.** rc7 (option 1), done; the rc6 package shipped rc6's form and disclosed the 11–20 µm.  RC8: `docs/ALIGNMENT.md` S1 quotes Xia's sentence (Review 1 of rc7, Q6): the text names rectangular-duct formulas (Yang et al.) without stating the shear coefficient, so it neither states nor rules out the parallel-plate form that panels d and e read.

### XIA-2. Panel e's registrations disagree across (V1, 3.7 px) (rc7 with XIA-1)

A consequence of XIA-1: the image's post columns on leg 3 are 19 µm closer than drawn, so the feature registration (the posts) and the chamfer registration (all edges) settle at different places across. Retest after XIA-1.  Done (rc7): on the parallel-plate mask panel e's registrations agree, V1 across +0.94 px (rc6 −3.73; `provenance/rc7_traces/README.md`).

### XIA-3. Panel f: leg 7's side channels read 20–33 µm narrower than drawn (open)

**What the figure shows.** A quarter of the way down leg 7, the side channels' outer walls lie inside the drawn ones: 23.6 / 19.7 µm on the printed reference, 29–33 µm on rc6's fig3 variant, 39.9 / 36.0 µm on rc7's aligned (parallel-plate) device, see INV-4 (`measurements/xia_side_walls_fg.json`, at panel f's own scale 1.774 µm/px; at the other panels' 1.798 the widths are 1.35 % larger, about 7 µm per side channel). At the outlet (panel g) the aligned device matches by construction (F).

**What has been tried.** F and D cannot match both panels (LIMITATIONS 24 d: lowering D fills the legs and widens panel f). The taper of XIA-1 widens panel f by another 7 µm (n = 1: 448 against 441 µm drawn there). So between panel f and the outlet the image's side channels widen faster than any constant-D schedule with Xia's layout.

**Possible causes (not settled by the figures).** Panel f's own scale (its post period is 1.5 % larger than the other panels'; if its true scale is the others', 7 µm of the 20–33 goes); Xia's loop pull-in applied to the next leg's recursion rather than to its widths (rc6's `loop_pull_mode: resistance`, the `_pullR` preset, is the nearest option; untested against panel f); a larger extraction on leg 7 than D implies.

**Status.** open: noted; retested after XIA-1 at rc7, with the pull-in mode as a candidate.

### XIA-4. Panel g's registrations disagree across by 1.84 px (note)

The outer wall and the posts read about 3 µm apart in relief contrast (V1 across 1.84 px, the gate 1.5 px). rc6 records it; no geometric change follows.

### XIA-5. Fig. 1e(ii)'s loop-1 widths: printed 69.5 / 57.6 / 77.4 µm, drawn 68.2 / 55.7 / 75.8 µm (note; rc7 with XIA-1)

The printed widths (leg 1's side channel before the loop; the loop's inner and outer side channels) are 1.3–2.0 µm wider than the printed reference draws. They are printed values (Fig. 1e caption), so full alignment asks for them. The aligned device and the XIA-1 taper do not close it (n = 1 draws 67.9 µm at leg 1's end). Xia's pull-in ("subtract 1–5 µm") and the loop's own width rule set the loop widths; rc6 takes 2 µm (**A**).

### XIA-6. rc6's search window against Xia's gap (rc7 window rule)

rc6's ± 8 px at 1.798 µm/px is ± 14.4 µm, nearly Xia's 16 µm gap (register GIF4-2). (Pre-rc7 numbers: the rc7 window is 7.19 µm, `provenance/rc7_traces/README.md`.)  At ± 3.5 px (± 6.3 µm) the aligned device's panels read medians of 1.2–2.8 µm and p90s of 3.5–4.8 µm, against p90s of 6.0–12.9 µm in rc6's window, where a post's edge can claim its neighbour's outline across the gap (`out/alignment.json`, the panel table). The posts' outlines lie about 2 µm outside the drawn ones along the flow (+1.9 to +2.5 µm per flank) and 2–4 µm across (rc6's V3): an outline bias of these fluorescence images in both directions. A narrower gap would move the flanks along the flow only, so the printed 16 µm gap stays. **Status.** rc7 window rule (with GIF4-2).

### XIA-7. Panels b, h, i are not traced (note)

Panel b is the inlet's serpentine head, which rc6 draws "as depicted" from it and from Fig. 3a (10 units measured, the rest extrapolated) but does not trace; panels h and i are the two collection ports. A trace of panel b (the head) belongs to rc7's trace code (not done in rc7); h and i show only port edges.

### XIA-8. Outside the model's validity envelope (note)

Both photomasks: u_min 0.022–0.023 below the κ plateau (the cutoff extrapolated, App. B C7 band), gap dominance M13 0.258 ≥ 0.05, and Re 60.6 > 10 at 800 µL/min (inertial flow at bends and gap mouths). rc6 reports all three; the geometry is unaffected (it follows f and the taper), the cutoff's absolute value is not calibrated there.

### XIA-9. Panels d (both masks) and g (the depicted mask): the outline metric covers under half of its samples (RC8, A7; recorded)

**What the traces show.** The outline metric's median and p90 describe the samples that found an image outline within
the window; rc7 did not count the others (the rc6 adjudication's A7).  RC8 counts them (`metric.coverage` = found /
eligible; the censored null ratios count a sample that found nothing at the window).  On Xia's relief-contrast panels
the walls are faint: panel d covers 0.42 (depicted) and 0.43 (aligned) of its eligible samples, panel g on
the depicted mask 0.48, and their censored null ratios are 1.00, 1.00 and 0.94: shifted by half a pitch or 5 µm across,
the mask scores as well, so the metric cannot tell the registration from a shifted one there.  (Viewed: on panel d the
overlay sits on the walls; the low coverage is the faint edges, which the metric cannot use.)

**Decision (D5).** Recorded, not claimed: the three panels' registrations (V1, V5) are reproduced and checked against
rc7's ranges, and their captions give the coverage.  Every claim now needs coverage ≥ 0.5 and is judged on the censored
null ratios where its limits use the metric nulls.  The chamfer nulls of panel d (the smaller chamfer null ratio 1.35 on the depicted mask, 1.37 on the aligned) could carry a registration claim stated as such;
panel g's across-flow chamfer null (1.05) cannot.

**Status.** RC8, done (`provenance/rc7_traces/code/make_anchors.py`, `tests/test_trace_anchors.py`).

---

## 5. Carried from the rc6 review adjudication (2026-10-18), § 12

All of these were planned as "rc6.1" or "rc7". With rc7 now the final freeze, they all go into rc7. The full wording is in `RC6_REVIEWS_ADJUDICATION_2026-10-18.md` § 12.

| id | what | kind | status at the rc7 freeze |
|---|---|---|---|
| A1–A3 | "Every panel …" overclaims: say which panels are traced (Xia 3b measured, not registered; 3h–i not traced; Dinh 3C left out) | doc | done 2026-10-19 (CHANGELOG, LIMITATIONS 24, the rc6 verification § 14) |
| A4 | the README sentence that runs the dfill and fig3 variants together; "inside" meaning narrower | doc | done 2026-10-19 (`provenance/rc6_traces/README.md`) |
| A4b, A5 | Gifford Fig. 4's along-flow tail: three texts disagree ("≈ 4 µm long", "7.5 µm on every post alike"); the data show a median 2.5 µm, mean +7.4 / +7.7, SD about 8 | doc | A4b done 2026-10-19 (verification § 14.2 notes); A5 moot: the README and LIMITATIONS texts it corrects were rewritten before rc7 (GIF4-2 holds the numbers) |
| A6 | the Dinh head corners, disclosed (superseded in substance by DINH-1's rc7 option; the disclosure stays for the arc default) | doc | done: LIMITATIONS 24 (f) for the arc default (2026-10-19); the depicted preset draws square turns (DINH-1); rc6's EVIDENCE.md generator not edited (rc6 history) |
| A7 | coverage: the metric's eligible and found counts in every record and table | record (trace code, not hashed) | partly: the rc7 traces README says the median and p90 are over the samples that found an outline (`metric.n`); the eligible count is not recorded (it needs the trace code changed and the traces re-run) |
| A8 | Xia Fig. 1e(ii): the printed loop-1 widths as an independent test (drawn 68.23 → 55.65 / 75.78 against printed 69.5 → 57.6 / 77.4 µm) | doc | done 2026-10-19 (PRESET_PROVENANCE loop pull-in; REPRODUCTIONS loop-1 row); rc6's EVIDENCE.md not edited |
| A9, A10, A15, A16 | provenance wording: Xia's scale arc, Gifford l_g's spread, "M" defined as image-measured, the glossary's words against the YAML keys | doc | done (A9 before rc7; A10, A15, A16 2026-10-19; A16 through `tools/options_sheet.py`) |
| A11 | EVALUATION: "validates every mask" → what is checked and what is not | doc | done 2026-10-19 (`docs/EVALUATION.md`; the tree-root copy refreshed from it) |
| A12 | the variant builds' commands in the traces README | doc | done 2026-10-19 (rc6 traces README; rc7's suite needs no variant builds) |
| A13 | figure titles and captions (clipped Gifford title; Xia caption on panels e and f) | doc | not done: rc6's composite figures; the rc7 packages show rc7's trace figures (`provenance/rc7_traces/results/`) |
| A14 | `panels_gifford6.py`: silence the all-NaN `nanmean` warnings | trace code (not hashed) | not done: a cosmetic warning in provenance code (results unchanged) |
| A17 | the verification's 68.3 µm against the plan's 68.23 | doc | done 2026-10-19 (the rc6 verification :576) |
| B1–B3 | rebuild the dfill variant record at the current stamp; ship the variant builds; a test that variant records carry the stamp | record | dropped: rc7 anchors no rc6 variant (`trace_anchors.json`); rc6's variant records stay history at their stamps |
| ADJ-C1 | `api.py` docstring: "the only physics inputs" → "the only targets" | code | done (rc7 code, 4a624626) |
| ADJ-C2 | the `rc6_values.json` note naming the arc value (hashed data) | code | done (rc7 code, 4a624626) |
| ADJ-C3 | the Dinh presets: "post pitch" → "row pitch"; rebuild the two Dinh records | preset + record | done (rc7 code, 4a624626) |
| ADJ-C4 | key aliases (`wall_along_flow_um`, `wall_cross_flow_um`) with a conflict check | code | done (rc7 code, 4a624626) |
| ADJ-C5 | a fluid section (ρ, μ, rheology) used by the gates | code | done (rc7 code, 4a624626) |
| ADJ-C6 | enumerated `assessment.status` values; `source_trace`; both in the CLI summary | code | done (rc7 code, 4a624626) |
| ADJ-C7 | test-suite hash, environment and lock hash in records | code | done (rc7 code, 4a624626) |
| ADJ-C8 | `defaults_used`: the effective value's origin | code | done (rc7 code, 4a624626) |
| ADJ-C9 | an optional passage report (straight sides, loops, transitions, necks against D and 2D), as a warning | code | done (rc7 code, 4a624626) |

## 6. Foundation, at the next unfreeze (not rc7)

D95–D97 (rc6) and D98–D105 (rc7: from the two foundation checks), as written in `foundation/ERRATA_NEXT_UNFREEZE.md`. rc7 does not change the frozen texts. If a device finding ever needed a frozen text changed, it would be listed here, not made in rc7.

## 7. The device packages themselves (not rc7)

| id | what | status |
|---|---|---|
| PKG-1 | PDFs carried a creation date, so reruns differed in three files | done 13:24 (no date) |
| PKG-2 | a run with the wrong Python stopped without logging why | done (package check first; STOPPED lines) |
| PKG-3 | a third figure: the whole device cut into a stacked sequence of rows (A. H. Routt, 13:27) | done (devicekit `render_stack`) |
| PKG-4 | DINH-P1 above | done |
| PKG-5 | a panel whose claim covers only some edges (Gifford Fig. 4: leg 4 claimed, leg 2 not) draws the unclaimed edges amber, as rc6's figures did | done 14:15 |
| PKG-6 | when `devicekit.py` is final, re-render every earlier package so all four carry the same kit (render step only; the figures already delivered are checked to stay byte-identical) | done 18:37: the final kit (`aa26bfa8…`) is in all four current packages (Dinh V1, Gifford Fig. 4 v2, Gifford Fig. 6 v2, Xia RS-7). Dinh re-rendered with it: every file in `out/` byte-identical except `run_log.txt` (and its line in `SHA256SUMS`); re-zipped. The Gifford version-1 packages are superseded by version 2 (PKG-9) and keep their kit |
| PKG-7 | the whole-device PNGs drew a rotated micrograph inside a black-filled bounding box (Gifford Fig. 4 at 45°: large black corners; Dinh at ≤ 0.6°: slivers, 2,304 pixels); the PNG path now gives the image an opaque alpha. PDFs unaffected | done 14:28 |
| PKG-8 | `registration.json` left the measured rotation empty for Gifford Fig. 6, whose records name it `rows_direction_deg`; devicekit now reads that key too (the same file in every package; the other packages' records use the other two keys, so their outputs do not change). Fig. 6 re-rendered: every figure byte-identical | done 16:59 |
| PKG-9 | Gifford Fig. 4 and Fig. 6, version 2 (A. H. Routt, 17:25–17:33: full alignment): the aligned build is each package's device and the rc6 design-level build its reference; a measure step checks the alignment (20 checks on Fig. 6, the table in each README); the package's first run writes the aligned photomask's sha256 and registration records into `reference/`, and a clean-room rerun from the zips checks them | done (Fig. 4: clean-room 18:25–18:36; Fig. 6: clean-room 18:36–19:09) |
| PKG-10 | Xia RS-7 package: rc6's fig3 variant as the device (it aligns Fig. 3a and panel g), the shipped preset as the printed reference; both photomasks and all twelve trace records byte-identical to rc6's; Fig. 3a drawn with the photomask on it | done (run 19:25–20:06; clean-room 20:40–21:30) |
| PKG-11 | `devicekit.py` additions for reproductions rc6 never recorded: `build_new` (a package reference the first run writes and every rerun checks), `check_build_physics` (network F to 10⁻⁴, no undeclared row outside the band, the declared departures exactly the expected ones, one component, nothing lost), `check_trace(establish=True)`; the trace code read from `rc6/trace_code` or rc6's own place `rc6/provenance/rc6_traces/code` (Xia) | done 17:46 (`aa26bfa8…`) |
| PKG-12 | the Xia package built its two photomasks one after the other in one Python process; the second differed from rc6's in three coordinates by 10⁻⁵ µm (BUILD-1). Its build step now runs each build in a Python process of its own, the two at once, as rc6 built them. The Gifford version-2 packages build their two in one process; each aligned photomask was checked to equal a build of its preset alone in a fresh process (Fig. 4 `3c72b2d7…`, Fig. 6 `7c77194f…`), and the design photomasks, built first, equal rc6's | done |

## 8. Found while building the packages (rc7 code)

### BUILD-1. A photomask depends on what the same Python process built before it (rc7)

**What happened.** The Xia package's first run built its two photomasks one after the other in one Python process. The first (the fig3 variant) equalled rc6's photomask byte for byte. The second (the shipped preset) did not: three of its 1,962,196 coordinates differed from rc6's by 1 × 10⁻⁵ µm, the last printed digit (sha256 `1249f737…` against rc6's `0161ab00…`; the same file size, code stamp and preset). Built alone in a fresh process, the same preset gives rc6's photomask byte for byte (19:07–19:24).

**Cause.** `engine/cifsched.py` memoises two quantities in module-level dictionaries:

- `bend_ratio`, the curved-duct conductance ratio, keyed by r1, r2 and h rounded to 10⁻⁴ µm;
- `depth_lambda`, the junction weight, keyed by h and w rounded to 10⁻³ µm.

Each cached value is the one computed at the first exact inputs that reached its key, so a later call with slightly different inputs receives the value for inputs up to half a rounding step away. At this device's sizes that is a bend ratio off by up to about 4 × 10⁻¹⁰ (relative) and a junction weight off by up to about 1.4 × 10⁻⁶ (absolute). Within one build the order of calls is fixed, so a build in a fresh process is deterministic: rc6's records, every package's clean-room run. A second build in the same process starts with the first build's entries. The two Xia presets share the depth and the loop radii, and three coordinates crossed a rounding boundary at the DXF's fifth decimal.

**Test (19:45–20:21).** The same two presets built in one process, in the same order, with only these two caches cleared between the builds (5,670 and 331 entries): the second photomask equals rc6's byte for byte (`0161ab00…`), the first equals rc6's too (`280a04ed…`). Without the clearing, the second build gave `1249f737…`. So the two caches are the whole mechanism.

**Physics.** None: 10 pm on three vertices of a 73 × 12 mm device, from resistances off by parts in 10⁸ or less. But a photomask should be a function of its input alone.

**Possible solutions.**

1. **(Preferred, rc7)** Clear both caches at the start of every `pipeline.run`. Every rc6 photomask stays byte-identical, because each was built in a fresh process, and a second build in one process then equals a build on its own. A test: two presets built in one process equal the same two built apart.
2. Key the caches by the exact inputs. This makes each value a pure function of its inputs, but changes rc6's photomasks in the last digit wherever two calls share a rounded key.
3. In the packages, build each photomask in a Python process of its own (done: PKG-12).

**Status.** done (rc7: option 1, with its test). Packages: done (PKG-12).

---

### BUILD-2. A leg spacing given too tight for a meander head drew the head across another leg; the build failed it only after drawing, and `check` passed it (rc7)

**What happened.** For the Fig. 4 folded-head variant (INV-2), a sealed 20 µm head was planned on the aligned preset's
layout: four legs in path order, 818 µm apart (the inset's spacing). The build drew the head's long first units across
leg 2 and the loop beside it, then failed the drawing on seven checks (part clearance 0.0 µm between "inlet and head" and
"leg 2", 50 walls lost at the floor, topology, side-channel continuity, the network read back from the drawing, its
cutoffs). A. H. Routt saw the crossing on the layout figure (2026-10-19 16:08: "This should not be possible … how can we
ensure that this never happens?").

**Cause.** The generator spaces a paperclip's legs to clear its head only when it chooses the spacing
(`layout._plan_paperclip`: max(2 hw_max, hw_head + hw_max) + 150 µm, here 1668–1672 µm). A spacing or leg axes that the spec gives
(`layout.leg_spacing_um`, `reproduction.leg_Y`) were used as given, and nothing compared them with the head before the
photomask was drawn and written:

- `opencif check` and `generate(check_only=True)` planned the design without laying it out, so they passed the spec;
- the collision checks ran on the written mask, so the failure came after a 14-minute build;
- the failing photomask was written under the deliverable name `<name>.dxf` (its record said FAIL).

A second, related defect: a spec value that a reproduction field replaces was only listed in the record's
`defaults_used` (code review R5b). So a `meander.head_width_um` given beside the preset's `reproduction.W` was not drawn,
and a retry with a narrower head drew the same head again.

**Fix (rc7, before its final freeze; A. H. Routt, 2026-10-19: "Fix in rc7 now").**

1. **Refuse colliding layouts before anything is written** (`build.collisions`, `build.plan_and_draw`). The build
   draws the outline in memory and checks it before the photomask is written:
   - every part (the inlet with the head, each leg, each U-turn, the outlet) at least the floor from every part it is
     not joined to;
   - the outlet ports separated;
   - every outlet port at least 2 mm from the inlet port.

   An outline that fails is REFUSED (exit 5, tag `layout`). The record names the two parts, the shortfall and the
   remedies (for a given paperclip spacing, the spacing the generator's own rule gives). When figures are drawn, a
   figure, `<name>_refused_layout.png`, outlines the two parts and their overlap. No photomask is written.
2. **`check` refuses what a build refuses.** `opencif check` and `generate(check_only=True)` draw the same outline in
   memory (`pipeline.layout_check`) and apply the same checks.
3. **The one self-placement collision the generator can cure, it cures; every other is refused.** Where only the
   outlet meets the inlet end, the outlet stretch is lengthened and the design re-planned (its closure takes the
   longer stretch). The case is a folded layout whose last leg ends beside its inlet with no loop there to clear (two
   nested legs: before, refused or failed; the fitter skipped such leg counts when it chose them), or a straight
   device too short for its ports. Not for a paperclip's ring outlet, or when another pair of parts collides too; it
   stops when a longer stretch does not help, and the gate refuses. With `legs: auto`, when every count collides, the
   fitter hands the build the first count it can cure. A design clear at once is drawn exactly as before. Any other
   collision of parts the generator placed is refused, never written.
4. **A failing photomask is never under the deliverable name.** The build writes `<name>_UNASSESSED.dxf`; its verdict
   names it: a FAIL (exit 2) `<name>_FAILED.dxf`, any other verdict `<name>.dxf` (`pipeline.name_mask`; arrays too).
   A run that stops before its verdict leaves no deliverable name. A composed array whose copies or manifold merge or
   disconnect is not built, and no photomask is written.
5. **Contradictory inputs are refused** (`spec.conflicts_with_reproduction`): a spec value that a reproduction field
   replaces with a different one, and a leg spacing given beside measured leg axes. The three Xia presets that
   depicted the head listed `meander.rectilinear_rows: 9`, which their depicted head (76 rows) replaced; the value is
   removed, and the head drawn is unchanged.
6. **A corner the sweep found** (`layout.leg_linear_local`): where a leg's last wall block runs past the leg's end
   into a U-turn, its corners' arcs could be cut at the leg's end, leaving a 0.19 µm step against the turn's wall that
   the corner-radius check FAILs (a short paperclip whose leg 1 holds only the head). The arcs now end on the band
   edges before the leg's end; only a leg whose last opening ends closer to its end than the corner radius changes
   (no shipped device: none ends a leg on a wall-block row).
7. **Evidence** (`tests/test_composition.py`, 22 cases):
   - 22 devices sampled from layouts, heads, inlets, outlets and start rules on four processes, at the plan stage,
     with full builds of a subset (8, of which 7 are drawn and 1 refused);
   - the generator's own placements are answered with their outline drawn (exit 0 or 3) and never refused for colliding
     parts; their outlet stretch is lengthened only where the outlet meets the inlet end; a spacing given too tight is
     refused with its remedies (a reconstruction of the case above; the case's own spec, at the final code, is refused
     by `opencif check` and `opencif build` in 50 and 66 s with no photomask written:
     `provenance/rc7_scripts/clearance_original_case/`);
   - every drawn build passes every hard geometric check (all but Q1 and Q2, the drawn mask's physics), with the exit
     code and the outlet stretch its check predicted; one case declares a Q2 FAIL (a floor start in evaluate mode, as
     m14);
   - a FAIL's photomask is renamed, and no `_UNASSESSED.dxf` remains.

An independent code review of the change (`review7/clearance/CODE_REVIEW4.md`, at code `0763a945e8f47c7b`: two major
findings, eight minor and two nits) found no way past the gate; its findings were fixed before the freeze (the conflict pairs it listed, the test's
build condition, the photomask's name before the verdict, the extension's limits and record, the check's spec, the
remedy wording), apart from the arrays' port spacing (F10), which LIMITATIONS item 29 states as not covered.  The
fixes were not re-reviewed as code; the final document check (`review7/final_check/FINAL_CHECK.md`, independent)
checked them against the diff and found two gaps, both disclosed: an array composed into a used folder keeps the
earlier photomask under the deliverable name (LIMITATIONS section E; `tools/make_arrays.py` now clears the folder, and
both shipped arrays were composed again), and two contradictory-input pairs the rule misses (LIMITATIONS 29).

**Effect on the shipped devices.** None: every one of the 50 built single devices (m14 included) passed the three
gated checks already (checked from their records). The rebuild at the final code `52ba91eddb3f07b8`
(rc7's `logs/compare_frozen_final.md`, in the rc7 evaluation archive): of the 53 records (51 builds, 2 arrays), 52 are the first freeze's, every number and
the photomask byte for byte where one is written (m14's under its `_FAILED` name); the 53rd changed:

- The Dinh figure-aligned refit (INV-1, D 3.67 µm) changes its own photomask, as intended.
(m14's deliberate FAIL keeps its photomask, byte-identical, under `_FAILED`.)

One record field changed on every build that runs the DRC: `validation.drc.file` (and the audit's `file`) names
`<name>_UNASSESSED.dxf`, the file the build audits before its verdict names it (item 4). Nothing reads the field; it was
found by the rebuild's comparison, after the freeze, and is disclosed in LIMITATIONS section E rather than fixed (a fix
would change the code stamp and need another full rebuild).

**Status.** done (rc7 final).

**RC8 (2026-10-21): the gaps rc7 disclosed, closed in the code** (the adjudication of the rc7 reviews, items K7–K11;
the RC8 code review).  K7: `spec.conflicts_with_reproduction` refuses the two pairs it missed, `meander.enabled` against
a contrary `reproduction.head` and an explicit head width (`narrow_um`, `channel_um`) beside a different
`reproduction.w_arm` (D3: a policy, one source for the head's widths).  K8: `pipeline.name_mask` names the file just
assessed and removes the other verdict's stale file, and `compose.array` clears its own outputs first, so an earlier
photomask can no longer stand under the deliverable name; every full build's record carries `validation.drc.file_final`.  K10: each
outlet-stretch extension try re-plans from the spec's own values, so a mask no longer depends on the tries before the
kept one.  K11: the 2 mm inlet–outlet rule applies between copies of an array (from each copy's inlet port as drawn, or
the shared port, to every other copy's outlet ports), gated before the array's photomask is written; the outlet ports
of neighbouring copies are reported, not gated (D4).  The full rebuild at code `593d977fc83ec217` reproduced every rc7
photomask byte for byte except the Fig. 4 reproduction's (`logs/compare_rc7_rc8.md`: of the 53 records (51 builds, one of them refused as in rc7, and the 2 arrays), 52 are rc7's in every number apart from what RC8 adds (`validation.drc.file_final` on the 50 full builds; the K13 note on `m14…_floorstart19`; `checks.inter_copy_port_spacing` on the 2 arrays), with 51 of the 52 photomasks byte-identical and none renamed; the one that changed is the Fig. 4 reproduction, with its input (D1: photomask `edaf92626d574aea…`, N 3543 as before)).  **Status.** done (RC8).

## 9. Investigations for rc7 (A. H. Routt, 2026-10-18 23:09: "review everything critically … regarding Dinh, rc7 should investigate it … consider anything else that might need further investigation")

### INV-1. Dinh's device length (DINH-3): the figures put the linear section at 514–642 rows, not 2,482 (rc7, decision)

Analysis `review7/inv1/`, checked independently (`review7/invcheck/INVESTIGATION_CHECK.md`, 2026-10-19 00:10: direction
confirmed, numbers corrected below).  The rc6 Dinh records and presets are unchanged.

**What the figures show.**

- **Fig. 2A is drawn to scale in both directions** (the inlet disc is round: b/a 1.011) **and its proportions are the device's.** Along the inlet port → middle outlet distance, the head starts at 8.9 %, the linear section runs from 54.2 % to 95.4 % (41.2 %); Fig. 2B's leader lines confirm the landmarks.  Fig. 3A's photograph puts the box C(i), "the beginning of the linear section", at 46–58 % (centre 52 %) of the same distance and C(ii) at 54 % of the linear section; perspective moves these by at most ~1.5 %.
- **Its length in mm depends on which printed footprint dimension the drawing spans** (its aspect is 18.6 : 1; the printed 3.2 × 75.8 mm is 23.7 : 1):
  - lateral scale (the arms span 3.2 mm; rc6's own `data/r4/dinh_layout_meas.json`): head 25.4 mm, **linear section 23.2 mm, 514 rows** at the 45.1 µm pitch;
  - length scale (the drawing spans 75.8 mm): head 31.8 mm, **linear section 28.9 mm, 642 rows**.
  - Whatever the scale, the photograph puts the linear section's start at ≥ 46 % of the chip's 75.8 mm, so the linear section is **at most ~36 mm (~800 rows)**.  rc6 draws 2,482 linear rows (112 mm) and a 25.7 mm head, 148.3 mm in all: 3–5 times the linear section the figures allow, and no 148 mm straight device fits a 75.8 mm chip.
- **The widths follow rc6's schedule shape, but cannot count rows.**  Fig. 3C's profiles read 222.7 → 245.8 µm wall to wall (profile scale; 227.7 → 240.8 at the pitch scale), against rc6's 225 → ~243.5 µm at the same fractions; the taper's width changes only ~10 % over the section, as a function of the fraction, whatever the number of rows.

**What the frozen physics says** (planner, `scripts/dinh_scan_D.py`; the check's `scripts/dinh_rows.py`).  [First analysis (00:10); superseded by the solution-1 numbers below (4.28–4.79 / 4.04–4.52 µm) and by codereview2 N4 (the head 24.1 mm from row 1, 1.3 mm short of Fig. 2A's 25.4 mm at the lateral scale; at the length scale Fig. 2A's head is 31.8 mm).]

- On the **lateral scale** the head at the printed f* (220 units, 25.7 mm in rc6's reconstruction) matches the layout's head (25.4 mm) almost exactly, and the linear section's 514 rows need a mean f per gap of **4.3–6.2 × f*** (F 0.505) or 5.0–7.2 × (F 0.559, the outlet identity at the printed widths).
- On the **length scale** the head is 24 % longer than f* predicts and the 642 rows need 3.7–4.9 × (F 0.505) or 4.3–5.7 × (F 0.559).
- Under the frozen law that is a linear-section cutoff of about **4–5.8 µm**, against 2.17 µm from the printed f* = 1.04 × 10⁻⁴.  So the figures fit **a head at the printed f* and a linear section at roughly 4–7 times f***, two cutoffs in one device (cleanest on the lateral scale, which rc6 itself recorded), or a head built differently from rc6's reconstruction.  A single cutoff fitted to the layout (≈ 4.5 µm everywhere) leaves the head about 3 times shorter than the layout's and is rejected.

**Performance (support, not proof).**  Dinh's filtrate (CIFside) carries PLTs of mean volume 6.8 fL (input 9.1 fL; Table 1) and 4.9 % of the 12 fL PLTs (p. 917).  As discs 0.66 µm thick (an assumption), 12 fL is ≈ 4.8 µm across and a 2.17 µm cutoff would pass mainly PLTs below ≈ 2.5 fL; a 4–5 µm cutoff passes them up to ≈ 8–13 fL.  The data favour the larger linear-section cutoff; platelets' effective size depends on their orientation.

**The meander and the length** (A. H. Routt's question, 00:23: "doesn't meander impact device length through resistance? Shouldn't meander be able to shorten device length?").

- **The row count is mass balance, not resistance.** App. A's march: P(i+1) = P(i)(1 − 2 f(i)), so 1 − F = Π(1 − 2 f_i). Given F and each row's fraction f_i (set by the cutoff through the lamina, f = g(δ*/w_c)), the number of gap rows is fixed; the side channels' resistance enters only through the f_i it makes each gap draw. Resistance decides what each side-channel segment must look like to make its gap draw f_i (the closure's R_s(k)).
- **What the meander does.** At the start the side flow S is small, so R_s(k) ≈ R_c P / S is large (exact under PE with equal segment lengths: σ = S/P = R_c/R_s, A:1148): a straight segment one row long would have to be narrower than the process allows. The head folds a long channel into each early row: a folded row needs Λ ≥ 2p (G:451; spec v5 G8, with the frozen p = GS + max(F_min, 2R_min) = 47 µm at the presets' 19 µm floor, so 2p = 94 µm; Dinh's 27 µm walls draw 2(GS + 27) = 110 µm of trunk per head row, against 45.1 µm per linear row). It makes a sealed start possible and lengthens each early row along the trunk; it cannot remove rows.
- **So the meander does set length, through the head.** The head's trunk length is at least N_fold × 2p, and N_fold, where folding stops, is set by the closure (spec v5 S:374, G11's handover): resistance sets length through how many rows must fold, not through N. At fixed N, D and F a tighter fold pitch, a finer handover width (Dinh: 19 µm, the floor, instead of the printed 21 µm folds 69 fewer rows, 4.5 mm shorter; foundation check 2) or a wider head budget (fewer arms, G5) shortens the device; nothing brings it below N × (G + g_p): the linear section is the shortest way to extract. An under-resistant fold over-extracts locally (a wider lamina, a cutoff above D) and the following rows under-extract by about as much; N and F do not move (spec v5 S:757; A:1218–1219). Only a design that runs the head above f* on purpose needs fewer head rows, by giving up D in them.
- **So the head's length measures the head's f.** The head ends at the seam, where a straight side channel at the row pitch reaches the printed w_S(in) = 21 µm (spec v5 S:374 with Dinh's printed width in place of the 19 µm floor: a G11 departure, disclosed); that needs a side flow σ_seam = r(w_c)/r(21 µm) = 0.0299 of the centre flow per side (w_c 75 µm, depth 140 µm; 0.0282 of the inlet flow per side, 5.65 % in all), fixed by printed widths, so the head holds ln(1 + 2σ_seam)/(−ln(1 − 2f)) rows: 279 at f* alone, 75 at 3.7 f*. rc6's head (18 measured rows at ~3.7 f*, which carry 0.0077, then 202 at f*) reaches it in 220 rows, 24–26 mm, the length of Fig. 2A's head; a head at ~4 f* throughout would be ~70 rows, 9–10 mm. Through resistance, the meander is what makes the head's drawn length a reading of its f.
- **The measured head agrees in part.** The 18 head units measured on Fig. 2A's inset (rc6 `data/depicted.json`) extract 2.8–13 f*, mostly about 3.7 f* (cutoffs about 4.2 µm), the scale the linear section needs; rc6 then fills the rest of the head's length with 202 rows at f*.

**Possible solutions.**

1. **[Superseded in its D by the revision below: D 3.67 µm, fitted to the printed length] (Built, rc7) A second, figure-aligned Dinh reproduction, `dinh2024_v1_figure_aligned`**, with no code change: the head exactly as the depicted preset plans it (prescribed as `reproduction.f_head`, 219 rows after code review R2: the 18 measured, then 201 at f*; `NM_force` 219) and the device's cutoff set for the linear section. The closure at the printed seam and outlet widths with the mass identity (A:1126–1132, A:1141, A:1148) fixes the linear section's total extraction (Σ −ln(1 − 2f) = 0.761 at F 0.5592), and App. A R13(a)'s invariant lamina (A:1294–1303, a premise: Dinh's rule, B:2958–2962) with the assumed taper turns Fig. 2A's 514–642 linear rows into 4.79–4.28 µm at the 1-D outlet identity's F 0.5592, 4.52–4.04 µm at the measured 0.505 (A:816–821). (App. A's inverse map, A:1158–1166, which an earlier version of this entry cited, needs every row's widths: corrected after foundation check 2.) Dinh's bead data (ESI Fig. S1B; App. B P-C7.2, B:2986–2996) put the cutoff below 4.7 µm, with soft edges (B:3069–3070); 96–98 % retention of the 4.7 µm beads favours the low end. The preset takes **D 4.29 µm (D, fitted** to Fig. 2A at the printed 75.8 mm length, 642 rows; it draws 641, 29.0 mm): a value inside the figures' range at either F (their overlap 4.28–4.52 µm) and the beads' bracket, at its low end. (rc7's first choice, 4.5 µm, the 1-D range's middle, was revised after foundation check 2.) The first build below was at 4.8 µm (514 rows), before the bead bracket was applied. Result (`review7/dinh2/`): exit 0, REPRODUCED with declared departures (seed; head drawn as depicted, rows 2–220, down to 0.45 D; the papers' transition); N 733 (513 linear rows, 23.2 mm), F 0.56006 (the outlet identity at the printed 44 / 52 µm), footprint 59.5 × 3.16 mm (Fig. 2A at the same scale: 60.6 mm; the chip: 75.8 mm). The linear rows extract 4.9 f* at 75 µm to 10 f* at 52 µm. Fig. 2B's four panels trace as the printed device's (they show local features both share); the discriminating evidence is Fig. 2A's and Fig. 3A's lengths. Its head rows at 2.17 µm lie below the κ plateau (the record notes it).  The shipped build, at 4.29 µm (code 4a624626): N 860 (641 linear), F 0.55964, 65.2 × 3.16 mm, exit 0.  Its two caveats (convergence review 1, M1 and M2): the head's agreement and the row count that sets D rest on different readings of Fig. 2A's scale (at the length scale Fig. 2A's head is 31.8 mm; its scale-free head : linear proportion is 1.09, the preset draws 0.83), and its 307 linear rows narrower than 62.7 µm (49 % of the section's extraction) run above κ_g's calibration cap 5.76 × 10⁻⁴ (f up to 8.36 × 10⁻⁴).
2. Keep the printed single-cutoff device (148.2 mm in rc7) and disclose: `dinh2024_v1_depicted`, the reference.
3. A single cutoff fitted to the layout: rejected (the head would be ~9 mm against 25 mm).

**What 1 contradicts.** Dinh's "to maintain the desired size cut-off throughout the device" (p. 916) and the printed f* in the linear section; App. B reads Dinh as one lamina throughout. It keeps every printed width, the depth, gap, GS, w_S(in), the outlet widths and the frozen physics row by row.

**Revision (rc7 final; A. H. Routt, 2026-10-19 15:45: "Fit printed 75.8 mm").** After the first freeze (4a624626) the two caveats were tested (`provenance/rc7_scripts/postfreeze/`):

- **The calibration (M2).** Gifford's own Fig. 5 devices show no drift of d50/D beyond the cap (5A 1.09–1.14, 5B 0.98–1.03, 5C 1.01–1.04, 5D 1.10–1.11), and the glossary's skim inversion holds to within 0.86 % up to f = 1.93 × 10⁻³ (G:1660–1663). A prefactor 10 % off beyond the cap moved the 4.29 µm reading to 4.13–4.51 µm. So M2 is a disclosed extrapolation, not a gap the foundation closes.
- **The beads** (`softedge.py`). On Gifford's measured retention curves (his Fig. 5, four devices of known D), the passage of Dinh's 4.7 µm beads (3.6, 6.7 and 14.2 % of the fluid's; ESI Fig. S1B) implies D 2.3–3.9 µm. At 4.29 µm Gifford's curves would pass 29–50 %; at 2.17 µm, 0–4 %.
- **The length** (`length_vs_D.py`). With the depicted head, the drawn device shortens by about 21 mm per µm of D near 3.7 µm. Dinh's printed 75.8 mm is reached at about 3.67 µm, moving about 0.09 µm per 2 mm.

**Decision.** The preset's D is fitted to the printed length: **D 3.67 µm** (`targets.cutoff_um: 3.67`, **D, fitted**). 4.29 µm stays documented as the reading from Fig. 2A's rows.

**Result** (the rebuild at the final code):
- N 1093 (874 linear rows, 39.5 mm), F 0.55927;
- 75.7 × 3.16 mm against the printed 75.8 mm;
- exit 0, REPRODUCED with declared departures;
- 88.7 % of F within 1 % of D;
- the beads' bracket holds: Gifford's curves pass 11–35 % at 3.67 µm.

**The caveats now.**
- **M1.** The printed length is kept at the price of the figures' proportion:
  - the linear section starts at 41.5 % of the inlet port → middle-outlet distance, against 54.2 % on Fig. 2A and 46–58 % on Fig. 3A's photograph;
  - it has 874 rows, against the photograph's bound of about 800;
  - the head : linear proportion is 0.65, against Fig. 2A's 1.09.

  No single reading of the figures gives the printed length, the depicted head and their proportion together.
- **M2.** 133 linear rows narrower than 53.6 µm run above the calibration cap. They carry 21 % of the section's extraction, with f up to 6.13 × 10⁻⁴ (`kappa_beyond_cap.py`: a prefactor 10 % off there moves D by −0.05 to +0.02 µm).

The fleet's bias statistic: on this preset the ±0.5 µm bias leaves 24 % of N outside D ± 1 % (was 29 %), but most of
those rows are its declared head rows, outside at nominal; the bias itself moves about 4 %, and the fleet's largest
bias share is 23 % (`dc_D3_F0p3`, `dc_D4_F0p3`).

**Status.** Both presets in rc7 (on A. H. Routt's instruction to continue until every goal is met; which is primary is his call); the figure-aligned one at 3.67 µm since rc7 final. The dissertation states the conflict either way.

### INV-2. Gifford Fig. 4: the inset's beads reject the aligned device's start (decision; RC8: the nested reading)

Analysis `review7/inv2/`, checked independently (`review7/invcheck/INVESTIGATION_CHECK.md`, 2026-10-19 00:10: conclusion
confirmed and strengthened; the remedy once proposed here, an open inlet, rejected by the same evidence; numbers
corrected below).

**The question.** The aligned Fig. 4 device (GIF4-1, decided 17:28; the v2 package's device and rc7's `gifford2014_fig4_aligned`) enters Gifford's schedule part-way: its side channels start at 51 µm, where the schedule has a side flow σ0 = 0.24 of the centre's per side, so about a third of the inlet flow (2σ0 / (1 + 2σ0) = 0.33) is in the side channels within the first rows.  App. A calls this start's σ0 "a fictitious inlet side flow".  Drawn sealed, the device pulls that fluid, with its beads, through its first 22 rows, whose cutoffs exceed 8.3 µm (rows 1–3 at 44–51 µm).

**Test 1, bead counts** (`scripts/particles.py`: a streamline model with steric exclusion on the aligned build's own row-by-row f and D; `beads_fig4.py`, `singles.py`): 33.2 % of the 8.3 µm beads should leave in rows 1–22 (about 29 % in rows whose cutoff is ≥ 20 µm, so it barely depends on the cutoff law), putting single large beads in the side channels at 0.42–0.45 (leg 2) and 0.30–0.34 (leg 4) of the centre's count.  The inset shows 1 and 1 (both, on the check's inspection, small-bead pairs or clusters; leg 2's side channels hold one clear large bead the classifier rejected; the classifier also undercounts the centre, which makes the test conservative).  Expected about 6 and 8.5 at 33 %; at 25 %, about 4.0 and 5.0.  From the counts alone the early loss is 6–8 % (maximum likelihood), 95 % upper limit 14–17 %.

**Test 2, the side channels' outer band** (the check's `scripts/lateral.py`, `outerband.py`; `out/outer_band_overlay_x3.png`).  In laminar flow each gap's fluid enters beside the posts and pushes older fluid outward, so the earliest-extracted fluid lies against the outer wall.  Under the 51 µm start, rows 1–3 send 15.2 % of the inlet suspension, at the inlet concentration of both bead sizes, into the side channels; at leg 4 it fills the outer 0.22 of the width (rows 1–22: the outer 0.39).  Measured on leg 4: the outer 0.2 of the width, holding the earliest 13.5 % of the inlet flow, should hold ≥ 98 beads at the inlet concentration and holds 2; the outer 0.3 should hold ≥ 150 and holds 15 (leg 2 the same way, more weakly: ≥ 26 vs 6).  This test needs no single/pair classification and does not depend on the 4.7 µm beads' behaviour near D.  **The earliest-extracted fluid was essentially bead-free** (Poisson p ≈ 10⁻³⁹ against the 51 µm start).  Combined estimate of the early loss: 0–3 %, 95 % upper ≈ 5 %.

**An open inlet is rejected too.** Side channels fed from the port would take S0 ≈ 24 % of the inlet suspension directly: count test p ≈ 2 × 10⁻³; outer 0.3 band ≥ 150 expected vs 15.  The image requires the side channels' first ~14–24 % of the inlet flow to be nearly bead-free: Gifford's design start (a near-zero first width at f: the first rows skim near-wall layers thinner than a bead radius, which the observed strong retention of 4.7 µm beads fits) or a side feed from a bead-free source, which the schematic does not show.

**Caveats.** Small counts in test 1; one snapshot (a start-up window of a few seconds is the only loophole the check found; settling and sticking do not explain the image); beads assumed spread evenly across the channel at the inlet.

**What it does to GIF4-1's options.** No reading fits every piece of evidence on Fig. 4:

| option | leg 2's widths | beads (no early loss; bead-free outer band) | schematic and caption ("second and fourth", nested) | single printed f_gap | App. A's start |
|---|---|---|---|---|---|
| A, Gifford's design start (rc6's `gifford2014_fig4_design`) | **20 µm narrow** per side | consistent | kept | kept | the sealed start (σ0 = 0) |
| C, start offset (the rc6 package's version 2; rc7 `provenance/rc7_traces/variants/gifford2014_fig4_start51.yaml`) | on (fitted) | **rejected** (33 % early loss; p ≈ 10⁻³⁹) | kept | kept | a fictitious σ0 |
| B, legs in path order (sequential; rc7 `gifford2014_fig4_aligned`) | on | consistent (rows 1–3 above D; 0.1 % of the large beads leave early) | **contradicted** | kept | the other start mode at Gifford's 10 : 1 floor (8 µm; 2σ0 = 0.0022); sealed only with a 4.86 µm first channel, 16.5 : 1 |
| D, a larger f over legs 1–2 (about +16 %, later rows about −15 %) | on | consistent (D about 5.6 then 4.8 µm, below 8.3) | kept | **contradicted** (two cutoffs) | sealed |

**Fig. 6 is unaffected.** Its start offset puts 10 % of the inlet flow in the side channels at the start (2σ0 = 0.11), consistent with the 80–85 % platelet recovery (Gifford's "yield", not App. B's F) (p. 4502); its micrograph carries no beads to test.

**Is there a reading with full alignment AND correct physics?** (A. H. Routt's question, 2026-10-19 00:23; built and traced on rc7, 00:44.)  No reading fits every piece of evidence (the table above); **reading B is the best-supported compromise**: every printed number kept, every visible feature of the inset aligned, the frozen physics exact on the drawing with a declared off-design entrance (corrected after foundation check 2, which found "Yes: B" and "strict physics" overstated).

- **What B is.** Gifford's printed design with every printed value (f_gap 3.34 × 10⁻⁴ at w_c 100 µm, depth 80 µm, final side width 300 µm, four legs) and his start rule (w_s(1) "close to or equal to zero, depending on the device production capabilities available", p. 4498) at the 10 : 1 aspect ratio he states (p. 4503): an 8 µm first side channel (**A**, an inference). The four legs lie side by side **in path order** (rc7 `layout.order: path`, 818 µm apart, **M**), so the leg beside leg 4 on the inset is leg 3. Leg 4's array ends 8.2 mm before the others (**D, fitted** to the neighbour's 124.7 µm through the row-fraction law; the schematic's inset box sits about a quarter of the way down the legs, the fit 0.19). G 17.2 / g_p 28.6 µm (**M**). Preset: rc7 `gifford2014_fig4_aligned` (C moved to `provenance/rc7_traces/variants/gifford2014_fig4_start51.yaml`).
- **Physics.** Every gap at D, except rows 1–3 (the 8 µm start; row 1 at 2.0 D, 0.5 % of the filtrate, declared): Gifford's start rule is App. A's "other start mode" (A:1118–1119; spec v5 S:271), a fictitious start flow 2σ0 = 0.0022 (230 times smaller than reading C's 0.485), realised on the sealed drawing as those three off-design rows. The 8 µm row-1 side channel is narrower than the 10.2 µm particles rows 1–3 pass. Exit 0, REPRODUCED with declared departures (floor start); every rc7 check passes; F 0.90650.
- **The image.** Traced with the rc7 code (window 8.37 µm, at most half the 17.2 µm gap): the claim's median 0.62 µm, p90 1.40 µm (n 7277; A: 1.09 / 2.51; C: 0.63 / 1.41); the registrations agree (V1 −0.09 / +0.03 px); the neighbour's four walls −0.77, −0.69, +1.04, −1.40 µm (A: +17.5 to +18.2); the legs' spacing 817.9 µm (818 drawn); the pitch 45.83 µm (45.8).
- **The beads.** The streamline model on B's own rows: only row 1 has a cutoff above 8.3 µm; 0.1 % of the large beads leave early, as the inset's nearly bead-free outer band requires.
- **The cost.** Three source items: the caption's "the second and fourth legs" (reading B: the third and fourth); the schematic's topology (its "not drawn to scale", p. 4502, is about scale, not which leg neighbours which; this register itself calls the topology unambiguous); and the same paper's Fig. 6, whose co-directional arrows show a nested layout. No printed number and no frozen statement is contradicted. (The schematic's inset box, a quarter of the way down its legs against the fit's 0.19, is a scale feature and is not counted as support.)
- **The alternatives, for the record.** Reading B sealed (the design-level start: every row at D, no declared rows): it aligns with the inset as well (the rc6 scratch trace of 2026-10-18 16:10: the neighbour's walls −0.78, −0.63, +1.09, −1.34 µm), since the start is not on the inset; but its 4.86 µm first side channel at 80 µm depth is 16.5 : 1 against the "10 : 1 aspect ratio limitation" Gifford prints (p. 4503), so it gives up a printed statement to avoid three declared rows (0.5 % of F): rc7 keeps the 8 µm start. A (nested, sealed): leg 2 20 µm narrow per side. C (51 µm start): rejected by the beads. D (a larger f over legs 1–2): contradicts the printed single f_gap.

**Decision taken (00:44, on A. H. Routt's instruction to continue "until complete, ensuring all goals are met"; reversible):** B is rc7's Fig. 4 reproduction; A stays the design-level reference; C is kept as a rejected variant in provenance. The dissertation states the caption/schematic conflict, as it states the paper's caption–text conflicts.

**RC8 (2026-10-21): the decision of 00:44 reversed by the taper direction (GIF4-1, RC8).**  Review 1 of rc7 measured
the inset's neighbouring leg's taper direction; opencif's own implementation (V13) confirms it: the neighbour flows the
way leg 4 does, so the legs are nested and reading B is rejected.  The table, with the taper column:

| option | the neighbour's widths | its taper direction (V13) | beads | caption, schematic, arrows | single printed f_gap | App. A's start |
|---|---|---|---|---|---|---|
| **nested, the 8 µm start, the inset's gap (RC8 `gifford2014_fig4_aligned`)** | **leg 2 17.5–18.1 µm narrow per side (open)** | **PASS** | consistent | kept | kept | the other start mode at 8 µm (2σ0 = 0.0022) |
| A, Gifford's design start (`gifford2014_fig4_design`) | leg 2 about 18 µm narrow per side | PASS | consistent | kept | kept | the sealed start (a 4.86 µm first channel, 16.5 : 1) |
| B, legs in path order (`gifford2014_fig4_pathorder`; rc7's device) | on (fitted) | **FAIL** | consistent | **contradicted** | kept | the 8 µm start |
| C, start offset (`gifford2014_fig4_start51`) | on (fitted) | PASS (nested) | **rejected** | kept | kept | a fictitious σ0 (2σ0 = 0.485) |
| D, a larger f over legs 1–2 | on | — | consistent | kept | **contradicted** | sealed |

So no reading fits every piece of evidence still, but the one that keeps every printed number, the caption, the
schematic and the measured direction is the nested one, at the price of leg 2's widths, which it states.  The beads'
analysis below is unchanged (it bears on the start, not the layout), and so is the 8 µm start.

**The start's width and the fabrication floor (rc7 final).** A. H. Routt's questions (14:38, 15:16): an 8 µm start lies well below a 19–20 µm floor; isn't that floor a limit of the photomask → mould → PDMS process? Can a simulation, sweep or validator close the gap? Decision, 15:51: "Keep 8 µm + variant".

- **The sources.**
  - Gifford 2014 printed his masks chrome-on-glass, adjusting the lateral design dimensions for the SU-8 feature size, and states "a 10 : 1 aspect ratio limitation" (p. 4503). Fig. 4's depth is "80 μm ± 3 μm" (caption) and "80 μm" (p. 4500), so 10 : 1 gives 8 µm.
  - Xia 2016: the original CIF devices' first side width was "limited by the minimum achievable feature size of the fabrication method", which the serpentine head removed. Its minimum feature is the gap, 16 µm, at 135–150 µm depth.
  - The Shevkoplyas lab's later guidelines (Lezzar et al. 2022, Sci. Rep. 12:13798): features "less than 20 µm are very difficult to mass produce", and an aspect ratio "above ~7:1 is very limited". At 80 µm these give max(20, 80/7) = 20 µm.
  - The lab's own devices pass both guidelines: Xia's 16 µm; Gifford Fig. 6's 17.1 µm gaps at 150 µm depth; Dinh's and Lezzar's 19 µm at 140 µm. **Fig. 4's own inset reads gaps of 17.2 µm (GIF4-3), below the 20 µm guideline**, so 20 µm was not this device's floor.
  - Gifford prints neither his first side width nor a feature floor: the 8 µm start stays an inference (**A**).
- **What the beads allow** (`provenance/rc7_scripts/postfreeze/early_loss.py`, a streamline model with steric exclusion; it reproduces this entry's 8 µm and 51 µm results). Each start was built on reading B's layout. The inset's outer band holds 2 beads; ≥ 98 would be expected at the inlet concentration, so the expected count is 98 × the early loss / 0.135.

  | first side channel | start flow 2σ0 | rows above D | early loss of 8.3 µm beads | beads expected in the outer band | P(≤ 2) |
  |---|---|---|---|---|---|
  | 8 µm (rc7, Gifford's 10 : 1) | 0.0022 | 1–3 | 0.087 % | 0.6 | 0.97 |
  | 11.4 µm (7 : 1) | 0.0076 | 1–6, 1.3 % of F | 0.65 % | 4.7 | 0.15 |
  | 12 µm | 0.0088 | 1–6, 1.4 % of F | 0.75 % | 5.4 | 0.095 |
  | 17.2 µm (the inset's gap) | 0.026 | 1–11, 3.6 % of F | 2.5 % | 18 | 2.8 × 10⁻⁶ |
  | 19 µm (the lab's later floor) | 0.035 | 1–12, 4.5 % of F | 3.4 % | 24 | 8 × 10⁻⁹ |
  | 51 µm (reading C) | 0.485 | 1–22 | 33 % | — | 10⁻³⁹ |

  (2σ0: the side flow the floor start assumes at gap 1, both sides, of the inlet flow, from each build's record.)  So the beads allow a straight start of up to about 12 µm and reject 17.2 µm and wider. A device held to a 20 µm floor can still start sealed: the rows narrower than 20 µm are folded into a meander head (spec v5's FOLD start, G11; Xia's remedy). The variant `provenance/rc7_traces/variants/gifford2014_fig4_fold20.yaml` draws it: a 59-row packing head, legs 1–2 at the generator's own spacing for the head (1672 µm; at the inset's 818 µm the head lay across leg 2, register BUILD-2), legs 2–4 on reading B's rows (`gifford2014_fig4_pathorder`, rc7's aligned preset; a path-order figure).  Built at rc7's final code, and again at RC8's (byte-identical, `b3fcda80…`): exit 0, N 3546, every check passing; no row's cutoff is above 8.3 µm, so the bead model loses none early.  Traced on the inset at its 20.5 µm gap: the claim's median 1.10 µm, p90 2.51 µm, as the design preset's (1.10 / 2.51 at the final code; reading B: 0.62 / 1.40 at the inset's 17.2 µm gap).  A figure, not a reproduction.
- **Status.** The 8 µm start is kept: the beads, Gifford's own stated limitation and the inset's 17.2 µm gaps support it. Disclosed in LIMITATIONS 28, PRESET_PROVENANCE and DISSERTATION_NOTES. The folded-head variant is shown as a figure, not as a reproduction.

### INV-3. Gifford Fig. 6: the platelet counts give F ≈ 0.69 in operation (note)

The typical experiment (p. 4502): input 384, centre output 1039, side output 93 (× 10³ µL⁻¹), yield 83.2 %. Mass balance: 93 F + 1039 (1 − F) = 384 gives F = 0.692 (side : centre 2.25 : 1); the yield equation gives the same. The caption's 140 µm gives 0.714 ("~2.5×"), the text's 144 µm 0.729. The operating F depends on the outlets' hydrostatic pressures (Gifford p. 4500: the ratio holds "provided the hydrostatic pressures in the output collection areas are kept equal"). The micrograph's width (145.5 µm at leg 7's end) remains the evidence for the drawn width (GIF6-2). Disclosed; no change.

### INV-4. Xia panel f (XIA-3): unchanged by the taper option (open)

Retested on rc7's aligned Xia device (`xia2016_rs7_aligned`, parallel-plate taper, D 6.4159 µm; `review7/xia1/`), each panel on its own registration (rc6's `xia_side_walls.py` method): panel f's side-channel outer walls lie **39.9 / 36.0 µm inside** the mask's (left / right; at panel f's own 1.774 µm/px, about 7 µm less at 1.798), against 33 / 29 µm for rc6's fig3 variant; panel g's wall −0.3 µm (fig3 variant −0.4). The taper option moves leg 3's centre onto the image (XIA-1) and fixes panel e's registration (V1 across 0.94 px, rc6 −3.73), but it does not explain leg 7's side channels: on leg 7 the centre is at its 150 µm floor in both forms, and D (6.42 against 6.57 µm) makes the side channels a few µm wider, not narrower. Traces at rc6's ± 8 px window: c 3.04 / 5.94, d 3.34 / 6.26, e 2.69 / 12.79, f 3.80 / 11.68, g 1.35 / 12.66 µm (median / p90; rc6 fig3: d 2.88 / 9.96, e 2.98 / 12.94, f 3.78 / 11.42).

**Status.** open, disclosed: the device's side channels on leg 7 a quarter of the way down read about 30–40 µm narrower (39.9 / 36.0 µm on the shipped aligned preset at panel f's 1.774 µm/px, about 33 / 29 at 1.798) than any constant-cutoff schedule that also matches panel g's last side channel. Tested since (rc7, `review7/inv4/`, planner only): Xia designing the side channels with the parallel-plate duct law too (r = 12 μ/(a b³), as the taper's parallel-plate shear might suggest). To draw panel g's 624.6 µm it needs F 0.893 (8.3 : 1) against the printed ~15 : 1, with D 5.87 µm to fill Fig. 3a's legs, and it still draws panel f's side channel at 436 µm (the exact law: 446.5; the image: 408–416): rejected. Panel f's placement by Fig. 3a's box is not the cause either (the image's 408–416 µm lies 179–226 rows, 9.3–11.7 mm, upstream on the drawn schedule, far outside the box). Candidates not tested: a side-channel law of Xia's own between panels f and g; panel f's scale (its own pitch). The rc7 package reports it.

### INV-5. Dinh Fig. 3C's scale bars (DINH-6): resolved (note)

The authors' intensity-profile axes in Fig. 3C confirm the 45.1 µm pitch scale (1.66 µm/px; the two panels' profile scales differ by 4 %, so each is good to about ± 2 %).  Fig. 3C's 100 µm bars, read with their anti-aliased ends (52.8 px), give 1.894 µm/px: **about 12–14 % short** (rc7's first reading, 16 %, took 52 px and Fig. 2B's bar at 46 px; with its ends Fig. 2B's bar is 46.95 px, 1.065 µm/px, and the pitch 45.04 µm, inside rc6's 45.1 ± 0.2).  Fig. 3C can now be traced at the pitch scale (rc7, trace).
