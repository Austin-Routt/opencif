# rc6 micrograph traces — the generator's masks on the papers' micrographs

opencif 1.0.0rc6, code stamp `4bbdaabef1225a24`; the masks are the shipped rc6 DXFs (`runs/`, sha256 in every record) and two
variant builds.  Nothing in this folder is hashed code.  The traces of rc5's masks (2026-10-16), from which rc6's
measured values came, are kept as records in `results/rc5_masks/`; their figures are in that delivery
(`OPENCIF_RC6_TRACES.zip` of 2026-10-16).

## What a trace is

An affine map from mask µm to image px, `img = origin + x·ex + y·ey`, never mirrored:

* **scale** — the panel's own scale bar *read with its anti-aliased end pixels counted by their coverage*
  (`trace_lib.bar_length_aa`): Dinh Fig. 2B 46.91 px / 50 µm, Gifford's Fig. 4 inset 107.52 px / 300 µm.  Where a bar
  cannot set a panel's scale, the panel's printed widths set it and the bar is a check: Xia's Fig. 3 panels b–g are at
  1.798 µm/px (`opencif/data/r6/rc6_values.json`: the mean of the printed centre widths, the post columns'
  centre-to-centre separations and panel d's loop curvature; the only bar, panel i's, gives 1.936), panel f at its own
  pitch (1.774 µm/px, its period 1.36 % off the others'); Gifford's Fig. 6 at 4.5248 µm/px from the printed 125 µm centre
  between the post rows' tip lines (`provenance/rc6_scripts/gifford_fig6_measurement.json`; its 250 µm bar runs off the
  raster and bounds the scale from above only, ≤ 4.5807).
* **rotation** — measured on the image (the post columns' slope; Gifford: the columns' principal direction or the rows'
  tip lines), never fitted to the mask.
* **translation** — fitted twice: on features (post centres matched to the mask's posts by a declared index rule from a
  non-periodic feature: the first post below a block, the last pill beside a wall, a printed box) and, independently, by
  chamfer matching of the mask's edges against the image's dark ridges.  Across the flow the two are independent; along
  a periodic array the chamfer scores one pitch away about as well as at the registration (V2's one-pitch null), so
  there the position is the non-periodic feature's (on Gifford's Fig. 6 verified row by row: END_CHECK).

## Validators (every panel, `results/*.json`)

* **V1** the two registrations agree (|Δ| ≤ 1.5 px along and across).
* **V2** nulls: the chamfer score and the nearest-outline metric of the same points shifted half a pitch along, one
  pitch along and 5 µm across.  In a periodic post array a one-pitch shift scores like the registration; the along-flow
  phase is fixed by the non-periodic feature the features were indexed from.
* **V3** signed nearest-outline offsets per class (post / block / island / channel × across / along × third of the
  panel); + = the image's outline outside the mask polygon (an island bigger, a channel wider).
* **V4** residual maps (`figures/panels/*_residuals.png`).
* **V5** scale scan: the mask rescaled about the region's centre by (s_along, s_across) ∈ 0.95–1.05 and re-registered;
  the pair the image prefers (1 = the declared scale).
* **V11, V12** (Xia's Fig. 3a, the whole-device drawing at its 5 mm bar): each part of the device registered on its
  own after the device is posed at its top (the shift of every group), and each traced panel's position against the
  box Fig. 3a draws for it.

The nearest-outline window is ± 8 px along each normal (8.5–22 µm on the panels at 1.07–2.79 µm/px); on Gifford's Fig. 6
(4.52 µm/px) it is ± 2.5 px (11.3 µm), so that it does not reach the next post across a 17 µm gap.

## Results (`results/`; the claims' limits in `trace_anchors.json`)

| panel | mask | µm/px | V1 along / across px | median / p90 µm | metric null ½ pitch / 5 µm | chamfer best / ½ pitch / 5 µm px | V5 along / across | role |
|---|---|---|---|---|---|---|---|---|
| Dinh 2B p1 head | shipped | 1.0660 | -1.09 / 0.18 | 1.29 / 2.70 | 1.75 / 1.82 | 1.53 / 3.23 / 2.37 | 0.980 / 0.988 | claim |
| Dinh 2B p2 transition | shipped | 1.0660 | 0.01 / -0.15 | 0.81 / 2.55 | 2.20 / 4.11 | 1.18 / 2.91 / 3.07 | 0.987 / 1.006 | claim |
| Dinh 2B p3 linear | shipped | 1.0660 | 0.00 / 0.11 | 0.44 / 0.89 | 2.71 / 3.84 | 0.57 / 2.61 / 2.48 | 1.000 / 0.995 | claim |
| Dinh 2B p4 outlet | shipped | 1.0660 | -0.18 / -0.09 | 0.56 / 1.57 | 1.82 / 4.20 | 0.95 / 2.44 / 3.16 | 0.996 / 1.001 | claim |
| Gifford Fig. 4 inset, leg 4 | shipped (design) | 2.7902 | -0.09 / 0.04 | 1.21 / 18.17 | 4.12 / 4.82 | 0.67 / 1.30 / 1.34 | 1.001 / 1.001 | claim |
| Gifford Fig. 6, leg 7 | shipped (design) | 4.5248 | -0.02 / -0.36 | 1.88 / 4.31 | 5.83 / 2.46 | 0.74 / 1.33 / 0.89 | 1.000 / 1.018 | claim |
| Gifford Fig. 6, leg 5 | shipped (design) | 4.5248 | 0.07 / -0.01 | 1.55 / 5.48 | 6.06 / 2.27 | 0.69 / 1.37 / 0.89 | 0.999 / 1.019 | claim |
| Gifford Fig. 6, leg 3 | shipped (design) | 4.5248 | 0.04 / 0.01 | 1.51 / 6.69 | 6.44 / 2.17 | 0.78 / 1.45 / 0.96 | 1.000 / 1.026 | claim |
| Gifford Fig. 6, leg 1 | shipped (design) | 4.5248 | -0.00 / 0.10 | 1.35 / 8.75 | 6.47 / 2.62 | 0.98 / 1.64 / 1.14 | 1.000 / 1.043 | claim |
| Gifford Fig. 6, legs 7-5-3-1 and their turns on one map | shipped (design) | 4.5248 | -0.01 / -0.10 | 1.47 / 5.93 | 5.99 / 3.06 | 0.72 / 1.36 / 1.00 | 0.999 / 1.001 | claim |
| Xia 3c leg 1, first posts | shipped (depicted) | 1.7980 | 0.13 / -0.01 | 3.06 / 5.96 | 5.49 / 2.95 | 1.26 / 3.01 / 1.74 | 1.000 / 1.003 | claim |
| Xia 3d legs 5, 3 below the loops | shipped (depicted) | 1.7980 | 0.26 / -1.24 | 2.28 / 9.05 | 8.10 / 4.81 | 1.84 / 3.54 / 2.47 | 1.005 / 1.009 | claim |
| Xia 3e legs 5, 3 half way | shipped (depicted) | 1.7980 | 0.22 / -3.67 | 3.03 / 12.93 | 10.69 / 3.77 | 1.76 / 3.90 / 1.96 | 1.000 / 1.008 | recorded |
| Xia 3f leg 7 (own scale) | shipped (depicted) | 1.7739 | -0.06 / 1.18 | 3.79 / 11.43 | 10.96 / 4.14 | 2.47 / 4.63 / 2.52 | 1.003 / 0.963 | recorded |
| Xia 3g outlet | shipped (depicted) | 1.7980 | 0.55 / 1.27 | 2.33 / 12.23 | 10.55 / 3.58 | 2.16 / 3.81 / 2.26 | 0.998 / 1.038 | claim |
| Xia 3c | variant dfill | 1.7980 | 0.13 / -0.01 | 3.05 / 5.96 | 5.48 / 2.94 | 1.25 / 3.01 / 1.74 | 1.000 / 1.003 | recorded |
| Xia 3d | variant dfill | 1.7980 | 0.19 / -1.03 | 2.34 / 9.41 | 7.90 / 5.00 | 1.86 / 3.52 / 2.47 | 1.006 / 1.008 | recorded |
| Xia 3e | variant dfill | 1.7980 | 0.22 / -4.20 | 3.58 / 12.80 | 10.67 / 3.76 | 1.87 / 3.95 / 1.96 | 1.000 / 1.009 | recorded |
| Xia 3f | variant dfill | 1.7739 | -0.06 / 1.18 | 4.09 / 13.00 | 11.28 / 4.56 | 2.47 / 4.64 / 2.52 | 1.003 / 0.973 | recorded |
| Xia 3g | variant dfill | 1.7980 | 0.55 / 1.27 | 2.30 / 12.16 | 10.55 / 3.58 | 2.16 / 3.81 / 2.26 | 0.998 / 1.038 | recorded |
| Xia 3c | variant fig3 (F, D off Fig. 3) | 1.7980 | 0.13 / -0.01 | 3.05 / 5.96 | 5.48 / 2.94 | 1.26 / 3.01 / 1.74 | 1.000 / 1.003 | recorded |
| Xia 3d | variant fig3 (F, D off Fig. 3) | 1.7980 | 0.28 / -0.71 | 2.88 / 9.96 | 7.45 / 5.14 | 1.99 / 3.54 / 2.49 | 1.006 / 1.008 | recorded |
| Xia 3e | variant fig3 (F, D off Fig. 3) | 1.7980 | 0.22 / -3.73 | 2.98 / 12.94 | 10.69 / 3.78 | 1.75 / 3.89 / 1.95 | 1.000 / 1.008 | recorded |
| Xia 3f | variant fig3 (F, D off Fig. 3) | 1.7739 | -0.06 / 1.16 | 3.78 / 11.42 | 10.96 / 4.12 | 2.46 / 4.62 / 2.51 | 1.004 / 0.950 (grid edge) | recorded |
| Xia 3g | variant fig3 (F, D off Fig. 3) | 1.7980 | 0.23 / 1.84 | 1.35 / 12.66 | 10.83 / 4.13 | 1.31 / 2.98 / 1.73 | 0.997 / 1.005 | recorded |

| whole device | pose (whole device) rotation / chamfer | far loops along / across | outlet ring along / across | V12 panels in their boxes |
|---|---|---|---|---|
| Xia Fig. 3a, shipped | -0.010° / 1.125 px | +5685 / +40 µm | +5405 / +36 µm | c ok, d ok, e ok, f ok, g OUT (8.2) |
| Xia Fig. 3a, variant dfill | -0.009° / 1.113 px | +94 / +41 µm | -186 / +36 µm | c ok, d ok, e ok, f ok, g OUT (1.1) |
| Xia Fig. 3a, variant fig3 | -0.008° / 1.111 px | +87 / +41 µm | -191 / +38 µm | c ok, d ok, e ok, f ok, g OUT (1.1) |

| Gifford Fig. 6 turns, mask | walls image / mask | image wall to nearest mask edge: median / p90 / max | turn centre mask - image, downstream / across | radii mask - image (about the image's centre) |
|---|---|---|---|---|
| uniform legs 588 µm apart (first form, variants/history588) | 16 / 18 | 23.6 / 59.7 / 230 µm | +11.5 / -1454 µm | — (not concentric with the image's turns) |
| the preset, with the measured leg axes (shipped) | 16 / 18 | 1.4 / 3.6 / 10 µm | +5.7 / +8 µm | +0.2 ± 2.6 µm |

The shipped preset's leg axes are read on these turns, so for it the nearest-edge distances, the centre across the rows and the radii agree by construction (consistency); the centre's downstream offset is the one independent turn result.  "(grid edge)": the V5 scan's best scale lies on the edge of its grid (a bound, not a minimum).

Roles: a claim passes its limits in `trace_anchors.json` (V1, the nulls — the metric's or, where relief contrast doubles the outlines, the chamfer's — the median and 90th percentile, V5); a recorded panel is reproduced and its ranges kept, not claimed.

**Dinh** (`dinh2024_v1_depicted`).  At the corrected bar the mask's pitch (G 19 + g_p 26.1 = 45.1 µm) is the image's
(45.08 / 45.09 / 45.16 µm on panels 4 / 3 / 2; V5 along 0.996–1.000 on the linear and outlet panels), and the outlet's
widths match the print and the mask to ≤ 0.5 µm (centre 52.44, printed 52; sides 44.08 / 43.94, printed 44; post
columns 106.06 µm apart, mask 106.0).  Panel 3 (row 972, 30 % along the linear section, located by the post columns'
centre-to-centre separation) reads w_c 68.8 and w_s 32.1 / 31.5 µm against the schedule's 68.1 and 32.1: the centre
agrees by construction (the row is the one whose post columns match the image's), the side channels are the test; the
outer walls would place the panel 87 rows earlier (row 885; they read 1.2 µm narrower than the mask's at row 972).  The
transition's blocks read about 2.2 µm shorter at each end along the flow than drawn (V3 block_along −2.10 / −2.47 /
−2.23 µm): the outline bias of a relief image, or blocks drawn 2 (s − GS) long where the device's are shorter — not
separable on one panel (LIMITATIONS 24 a).  The head panel (p1) prefers a 2 % smaller scale along and across (V5 0.98 /
0.99); its units are the depicted head's (`opencif/data/depicted.json`).

**Gifford, Fig. 4 inset** (`gifford2014_fig4_design`, leg 4's end).  Pitch 45.83 / 45.84 µm (mask 45.8); post columns
159.2 µm apart (mask 159.5 = 100 + 59.5); median 1.21 µm; V5 1.001 / 1.001.  The neighbouring leg (leg 2) lies 817.9 µm
away on the image.  Gifford prints no spacing, and the generator's equal-leg rule gave 969 µm; rc6's presets take the
inset's (818 µm; `provenance/rc6_scripts/apply_rc6_layout_presets.py`), so leg 2's centre line and posts sit on the
image's -- by construction, consistency rather than a test.  Leg 2 stays in the no-claim colour: its side channels there
read 17–18 µm wider per side than the rows the mask places at that point (`measured.leg2_side_walls`), and which of leg 2's
rows pass the inset is set by the legs' lengths, which the inset does not show.  The claim's 90th percentile, 18.2 µm,
is the posts' long sides: the image's outlines lie about 7.5 µm outside them along the flow (V3 post_along +7.4 / +7.7
µm, n 356 / 1344), on every claimed post alike, while the posts' centres coincide (V1 −0.09 / 0.04 px; the feature
fit) -- the metric's outline detection on this relief image (the dark band beside each pill's flank), not a geometric
offset (the
metric's nulls carry the same p90, 18.2–20.8 µm); its limit (20 µm) was set at this result.

**Gifford, Fig. 6** (`gifford2014_fig6_design`; new in rc6).  The micrograph shows the downstream ends of four legs and
the start of three turns.  *Which legs* (`panels_gifford6.py`): the caption's arrows give one flow direction for all
four (towards the walls; read on the figure), and in an opencif paperclip odd legs flow +x, so the windows of four
neighbouring legs left are (7, 5, 3, 1) and its mirror; computed: the side channels at the four ends, read on the image
without the mask (wall line to the posts' outer tip line), are 145.5 / 114.1–114.5 / 87.1–87.5 / 58.5–59.3 µm against the
drawn 139.9 / 108.7 / 80.7 / 49.3 µm, an offset of one sign on every leg (+ 7.0 ± 1.7 µm; per leg 5.5 / 5.6 / 6.6 / 9.6 µm, growing as the
channel narrows) for (7, 5, 3, 1) and none for the mirror (SD 60 µm) -- of the windows the flow direction excludes, the
mirrored (6, 4, 2, 1) also shows a near-uniform offset (17.5 ± 5.1 µm), so the flow direction read on the figure is what
excludes it; and, seen on the figure, the top leg's channels run straight on (the outlet) while the others curve (the turns).
*The layout* (rc6, after the author's review of 2026-10-17).  Gifford prints none (the schematic is not to scale, p.
4502).  The generator's uniform paperclip (723 µm; 588 µm, the image's leg spacing, in a first variant kept in
`variants/history588/`) turns each leg about a point one half-spacing beyond leg 1, far tighter than the image: its 16
turn walls of legs 5, 3 and 1 are arcs about one common centre 1745.7 µm beyond leg 1 (residual SD 0.13 px, 0.6 µm;
the centre's SD 2.2 µm in the common fit, but 26 µm between the centres the walls give one by one, and legs 2,
4, 6 at 2 Y_c − Y carry twice the centre's error; radii 1.68–3.17 mm; `code/gifford6_turns.py`, `measurements/gifford_fig6_turns.json`,
`figures/fig_gifford2014_fig6_turns.png`), each turn carrying the centre channel (126 µm) and both side channels as
concentric bands, as opencif draws them.  The presets therefore take the image's layout (`reproduction.leg_Y`,
`reproduction.pc_last_short`): legs 1, 3, 5, 7 on the post rows' tip lines (0 / −587.8 / −1175.4 / −1764.3 µm), legs 2,
4, 6 at 2 × 1745.7 µm − Y of the leg before (3491.4 / 4079.3 / 4666.8 µm) -- semicircular turns about the measured
centre, opencif's turn model; only 8–23° of each turn is in the frame, so that closure, and with it legs 2, 4, 6, is an
assumption -- and the last leg's array 257 µm short.  *Traced:* each leg on its own map (V1 ≤ 0.36 px; medians
1.4–1.9 µm), and the four legs with the turns of 5, 3 and 1 on ONE map fitted on the ends of legs 5, 3 and 1 (V1
0.01 / 0.10 px; median 1.47 µm, p90 5.93 µm; V5 0.999 / 1.001).  The turns compared directly
(`panels_gifford6.turn_compare`: every point of the 16 walls followed on the image against the mask's nearest edge,
no window; consistency for the turns' shape, since the leg axes, hence the radii and the centre across the rows,
were read on these walls): median 1.4 µm, 90th percentile 3.6, max 9.5 (the uniform 588 µm layout: 23.6 / 59.7 / 230 µm, its
turn centre -1454 µm across from the image's); about the image's turn centre the mask's walls lie at radii +0.2 ± 2.6
µm from the image's (the side channels' own offset (a) included); fitted alone, the mask's turn centre lies +8 µm
across from the image's.  The one independent turn result is along the flow, which nothing read on the turns sets
(the generator places the turns from the legs' ends): the mask's turn centre lies 5.7 µm downstream of the image's
(the image centre's SD there 0.5 µm).  Leg 7, not in the fit, lies 0.31 µm across from where its own end registers it --
consistency, since the axes were read on this image (`figures/fig_rc6_fix_gifford6_turns.png`: before and after).
*What is tested and what is not:* the scale is set by the printed centre, and the pitch, the post length and the leg axes
of these presets were measured on this micrograph, so the centre, the pitch (V5 along), the posts' across positions and
the legs' positions agree by construction; along the flow a periodic array scores one pitch away as well as at the
registration (chamfer and metric), so the along position rests on pairing each row's last pill with the mask's last
post, which **END_CHECK** verifies on both sides, row by row: on the image the last pill is followed one gap downstream
(16.7–17.9 µm) by the divider wall's end and by no further outline for a pitch; on the mask the last post is followed
20.5 µm downstream (G) by the divider wall.  Tested, then: the leg identity, the leg ends, the side channels and the
turns' start along the flow; the turns' shape, given axes read on them, is consistency.  Two things differ, each
reported: (a) the side channels at the leg ends are 5.4–10.0 µm wider on the image (V3 channel_across +4.5 to +9.0 µm per wall; V5 across 1.018–1.043 on the per-leg maps): the
assumed gap (G 20.5, Fig. 4's; this image reads 17.1 µm along the flow), the outline bias of straight walls against
rounded tips, or the device itself -- not separable here; the final side channel reads 145.5 µm (caption 140, text 144),
but an offset of one sign on every leg makes that a systematic reading, not a choice between the two printed values
(the offset is not constant, growing from leg 7's 5.5 µm to leg 1's 9.6 µm as the channel narrows, so the image does
not fix the final width to better than a few µm) (the metric's
± 11.3 µm window drops part of leg 1's walls, which lie ~9 µm out); (b) leg 7's drawn end is 4 rows beyond the image's,
because each middle leg is drawn a whole number of rows and the last leg takes the remainder (LIMITATIONS 25).

**Xia, Fig. 3** (`xia2016_rs7_depicted`, at 1.798 µm/px).  Panel c: leg 1's centre at its first posts 299.9 µm (printed
300).  Panel d: legs 5 and 3 leaving the top loops, 150.3 and 219.8 µm against the mask's 150.0 and 232.7, their axes
902 µm apart (mask 905).  Panel e: half way down legs 5 and 3, at Fig. 3a's box: 151.5 and 186.0 µm against 150.0 and
205.1 — leg 3 is narrower than the constant-shear taper puts it, so one translation cannot fit both legs (V1 across
−3.7 px: recorded, not claimed).  Panel f: leg 7, 151.5 µm (mask 150.0).  Panel g (the outlet): V5 across 1.038, the
last side channel 624.6 µm against the schedule's 600.  The taper's form is **A** (Xia prints only its ends): the
micrographs show the device narrowing faster than the constant-shear assumption over legs 2–4 (LIMITATIONS 24 c).
On panel g the chamfer's 5 µm-across null barely discriminates (2.26 against the best 2.16 px; its limit, 1.04 ×, was
set at this result): there the metric's across null (3.58 against 2.33 µm) and V1 carry the claim across the flow.
The relief contrast doubles many outlines, so the metric's 5 µm null is weak on these panels; the chamfer's nulls carry
the claims.

**Xia, Fig. 3a — the whole device** (`xia_fig3a_device.py`).  At its 5 mm bar (178.32 px) and a pose fitted at the
device's top, the shipped mask's legs run across within 2–41 µm of the drawing's, but the drawing's far loops and its
outlet ring lie 5.69 and 5.41 mm beyond the mask's (V11): at the pitch the micrographs show, Fig. 3a's legs hold 1190
rows where the schedule places 1079.  V12: panels c–f lie inside their boxes; g lies 4.9 mm from its box (the same
deficit).  *The variant* `variants/xia2016_rs7_depicted_dfill.yaml` (D 6.5315 µm, × 0.9515; F unchanged) fills the
legs -- with N pinned by the drawing's length, two targets meet one free variable (App. B, C8 CQ5, `pinned_row_count`),
so this variant keeps F and lets D go (the fig3 variant below reads both off the figures).  Registered the same way, the
dfill mask's far loops lie 0.09 mm and its outlet ring 0.19 mm from the drawing's, and its legs run across within
2–41 µm as before.  V12 (positions along Fig. 3a, not widths): panels c–f lie inside their boxes; panel g lies 0.65 mm
from its box's centre (1.08 of the box's half-size), 0.46 mm of which is panel g's own offset from its box, which the
shipped mask shows as well (a 5.41 mm ring shift against a 4.94 mm panel offset).  Its panels trace as the shipped
mask's (medians 2.3–4.1 µm).  Filling the legs by F
instead (0.958) would end the last side channel at 860 µm against the 624.6 µm panel g shows, so the micrographs choose
the cutoff, a conditional like Dinh's length (LIMITATIONS 24 d).  *The variant* `variants/xia2016_rs7_depicted_fig3.yaml` reads both decoupled targets off Fig. 3: F 0.94072 (15.9 : 1), at which the schedule's last side channel is panel g's (624.6 µm drawn, 624.6 µm on the image; the printed ~15 : 1 draws 600), and D 6.5705 µm (× 0.9572 of the printed f_ref's), at which its 8106 rows fill Fig. 3a's legs at that F (`code/xia_fig3_fit.py`: F for the width, D for the rows, alternately; `results/xia/xia_fig3_fit.json`).  Panel g: V5 across 1.005 (the preset's 1.038), median 1.35 µm (the preset's 2.33); its two registrations now differ by 1.84 px across (the preset's 1.27; V1's limit 1.5): with the side channel matched, the left wall now pulls the chamfer, and the wall and the posts disagree across by about 3 µm, the size of the outline bias of relief contrast (the posts read +2.3 / +0.8 µm against the drawn width, the centre 1.7 µm narrower; the 624.6 µm target is outline to outline, and a 1-3 µm correction would raise F by 0.0001-0.0003); recorded, not claimed.  On Fig. 3a its far loops lie +0.09 mm and its outlet ring -0.19 mm from the drawing's; V12: c inside, d inside, e inside, f inside, g 1.08 of the half-size; its panels trace as the preset's (medians 1.3–3.8 µm).  `figures/fig_rc6_fix_xia3g_wall.png`: panel g before and after.  Panel f (leg 7, a quarter of the way down; recorded, not claimed) reads the other way (`code/xia_side_walls.py`, `measurements/xia_side_walls_fg.json`, each panel on its own registration): its side-channel walls lie 24 / 20 µm nearer the centre line than the preset's (left / right) and 33 / 29 µm nearer than the variant's (the image's side channels are narrower than both masks draw; rc7's parallel-plate aligned preset: 39.9 / 36.0 µm, register INV-4) at panel f's own 1.774 µm/px (from its pitch) -- about 7 µm less at the 1.798 of the other panels -- and the variant's scale scan across on panel f ends on its grid's edge (0.950); at panel g the image's wall lies 26 µm outside the preset's and 0.4 µm from the variant's.  If panel f's scale holds, the device's side channels widen faster along leg 7 than the model's and no F or D matches both panels (F sets the outlet's width, the outlet identity): the schedule's form would depart from the assumed constant-shear one, as the centre taper's does (LIMITATIONS 24 c); one panel at its own scale does not settle it.  Both variants are readings of the figures under the model's premises, not of the text; the shipped presets keep the printed values.

## Whole-device figures (`full_device/`)

`full_device.py` draws every polygon of a shipped DXF with the micrographs placed by the inverse of the maps their traces
stored (nothing re-fitted; V9 checks the round trip, V10 the single scale), and, for Xia, Fig. 3a laid under the device
by its whole-device pose: `<name>_full.pdf` (vector), `<name>_strip<i>.png` (2 µm/px), `<name>_zoom_<letter>.png` and
`full_device_composite.png` (the four shipped devices at one scale, 4 µm/px: their strips reduced by 2 with a box
average, x0 and y0 kept on pixel edges, `code/composite_from_strips.py`; drawing all four in one figure needs more
memory than the container has), with `full_device_index.json`.  Drawn: `dinh2024_v1_depicted` (148.3 × 3.2 mm), `gifford2014_fig4_design` (47.4 × 5.0 mm), `gifford2014_fig6_design` (53.8 × 8.3 mm), `xia2016_rs7_depicted` (73.3 × 12.1 mm), `xia2016_rs7_depicted_dfill` (78.9 × 12.1 mm), `xia2016_rs7_depicted_fig3` (78.9 × 12.2 mm); V9 and V10 hold for every view (0 failures); the composite holds the four shipped devices at 4 µm/px.  To regenerate the folder for any
tree: `python3 -I code/run_full_device.py <repository> <rasters> <out> [<variant builds>]` (the configuration, one
process per device, the merged index, the composite); the shipped `config*.json` and the figure specs
(`figures/*.json`) are this run's records and carry the build machine's paths, which those scripts rewrite for yours.

## Files

* `code/` — `trace_lib.py`, `panel_core.py` (one panel: chamfer, V1–V5, figures), `panels_dinh.py` + `panels_dinh_v2.py`
  (Fig. 2B p1–p4), `panels_gifford.py` (Fig. 4), `panels_gifford6.py` (Fig. 6), `panels_xia.py` (3g), `panels_xia_more.py`
  (3c–f), `xia_fig3a_device.py` (Fig. 3a: V11, V12), `xia_arc_scale.py`, `xia_leg_fill.py` (the D that fills Fig. 3a's
  legs), `xia_fig3_fit.py` (F from panel g's last side channel and D from Fig. 3a's legs together), `gifford6_turns.py`
  (Fig. 6's turn walls and their common centre), `xia_side_walls.py` (panels f and g's side-channel walls against each mask), `alignment_fixes.py` (the three before / after figures), `measure_scales.py`, `full_device.py` + `full_device_config.py` + `composite_from_strips.py` + `run_full_device.py` (the three in order), `composite.py` + `composites_rc6.py` (the
  dissertation figures redrawn from the records), `summary_table.py` (the table above), `evidence_index.py` (EVIDENCE.md), `make_anchors.py`, `run_*.py`
  (`run_gifford6.py` for Fig. 6), `build_variant.py`, `zoom.py`, `dinh_length_at_corrected_pitch.py` and
  `apply_doc_edits.py` (rc5-era).
* `results/` — one JSON per panel and mask (maps, V1–V5, widths, captions); `results/rc5_masks/` the 2026-10-16 traces.
* `figures/` — `fig_dinh2024_fig2B_traced_rc6`, `fig_gifford2014_fig4_inset_traced_rc6`, `fig_xia2016_fig3_traced_rc6`
  (PDF, PNG, and the composite spec JSON), `xia_fig3a_device_{rc6,dfill,fig3}.png` (the whole-device registrations),
  `fig_gifford2014_fig6_turns.png` (the turn walls and their common-centre circles), `fig_rc6_fix_gifford6_turns`,
  `fig_rc6_fix_xia3g_wall`, `fig_rc6_fix_gifford4_start` (before / after, `alignment_fixes.py`);
  `panels/` each panel at 600 dpi with its caption (PNG, PDF), `_clean.png`, `_residuals.png` (Gifford's Fig. 6 also
  `gifford_fig6_legs__rc6`, the four legs each on its own map).
* `measurements/` — `scales.json`, `xia_fig3d_arc_scale.json` (panel d's loop arcs paired with the rc6 mask's walls: 1.792 ± 0.004 µm/px; the rc6 scale's curvature term uses only loop 4's axis radius against Fig. 3a's axes, 1.794, independent of the mask), `gifford_fig6_turns.json`, `xia_side_walls_fg.json`, `summary_table.md`.
* `variants/` — `xia2016_rs7_depicted_dfill.yaml`, `xia2016_rs7_depicted_fig3.yaml` and their records (`records/`);
  `history588/` the first form of rc6's Fig. 6 layout (`gifford2014_fig6_design_fig6legs.yaml`: uniform legs 588 µm
  apart; superseded by the presets' measured leg axes), kept with its record for the before / after figure; the DXFs are not in the repository (`python3 provenance/rc6_traces/code/build_variant.py .
  provenance/rc6_traces/variants/<name>.yaml provenance/rc6_traces/variants/runs/<name>`, ≈ 30 min each;
  `history588/` likewise, ≈ 10 min); `variants/rc5/` the two rc5-era variants (`gp26`, `s181`), whose values rc6's presets now carry.
* `EVIDENCE.md` — the evidence figure by figure: preset, zoomed trace, whole device, numbers, and what is tested and
  what agrees by construction (`code/evidence_index.py`, from the records).
* `trace_anchors.json` — inputs, scales, each panel's anchored values and limits, Fig. 3a's registrations, the
  variants' declared fields: `tests/test_trace_anchors.py`.

## Reproducing

The rasters are the publishers' and are not in the repository.  The JSON configurations in `full_device/` and
`figures/` name their inputs by paths relative to the rc6 working folder (`rasters/` for the page images below,
`opencif_rc6/runs/` and `variants_out/` for the rc6 builds); every raster and photomask is identified by its sha256 in
the records, which the scripts check.  From the three PDFs:

```
pdfimages -png -f 5 -l 5 "Dinh 2024 ….pdf" dinh5     # dinh5-000.png  sha256 a21e2e3ce0710e32…
pdfimages -png -f 5 -l 5 "Gifford 2014 ….pdf" gif5   # gif5-001.png   sha256 e70a8b38b87f9639…
pdfimages -png -f 7 -l 7 "Gifford 2014 ….pdf" gif7   # gif7-000.png   sha256 182826e7aaccf14c…
pdfimages -png -f 5 -l 5 "Xia 2016 ….pdf" xia5       # xia5-000.png   sha256 e9fb7fbb348567d3…
export OPENCIF_PAPER_RASTERS=<that folder>
python3 -I provenance/rc6_traces/code/build_variant.py . provenance/rc6_traces/variants/xia2016_rs7_depicted_dfill.yaml <builds>/xia2016_rs7_depicted_dfill
python3 -I provenance/rc6_traces/code/build_variant.py . provenance/rc6_traces/variants/xia2016_rs7_depicted_fig3.yaml <builds>/xia2016_rs7_depicted_fig3
export OPENCIF_TRACE_VARIANTS=<builds>
python3 tests/test_trace_anchors.py            # 188 checks, 0 failed (exit 0); without the variant builds: 134 checks, 0 failed, 2 inputs missing (exit 4)
# (rc6's suite.  In the rc7 tree the same command runs rc7's anchors, 268 checks, without variant builds: provenance/rc7_traces/README.md)
```

On rc6 (code `4bbdaabef1225a24`, 2026-10-17, after the review of the layout work): 188 checks, 0 failed, 0 inputs missing, 43 min on two cores shared with the batch.  The acceptance limits were set with the results in view: they are regression limits for rc6, not a
blind test.  What is blind is each panel's scale, its rotation, and the agreement of two independent registrations.

To regenerate the records and figures instead of checking them (each script prints its usage; ROOT the repository,
RASTERS the folder above, OUT a results folder, TAG `rc6` for the shipped masks, and `PRESET DXF RECORD LABEL` after
the TAG for any other mask, such as a variant built by `code/build_variant.py`): the panels with `code/run_dinh.py
ROOT RASTERS OUT TAG`, `code/run_gifford.py ROOT RASTERS OUT TAG`, `code/run_gifford6.py ROOT RASTERS OUT TAG
[--one-map --turns]`, `code/run_xia.py ROOT RASTERS OUT TAG rc6` and `code/run_xia_more.py ROOT RASTERS OUT TAG cdef`;
the measurements with `code/xia_fig3a_device.py`, `code/gifford6_turns.py`, `code/xia_arc_scale.py` and
`code/xia_side_walls.py`; then `code/summary_table.py`, `code/make_anchors.py`, `code/composites_rc6.py`,
`code/alignment_fixes.py`, `code/run_full_device.py` and `code/evidence_index.py`.  The masks themselves come from
the package (`python3 -m opencif build presets/<name>.yaml`, or `tools/run_batch.py` for all of them).

Licences: Xia 2016 is CC BY 4.0 ("This work is licensed under a Creative Commons Attribution 4.0 International
License.": attribute, and say the figure is adapted).  Dinh 2024 and Gifford 2014 carry "© The Royal Society of
Chemistry" with no Creative Commons statement in the PDFs: RSC asks for a request through the Copyright Clearance
Center (free for a thesis submitted to the author's institution: choose "Academic institution") and the line
"Reproduced from Ref. XX with permission from the Royal Society of Chemistry" (for an overlay, "Adapted from …")
(rsc.org, Licences, copyright & permissions, read 2026-10-16).  The rasters are not redistributed here.
