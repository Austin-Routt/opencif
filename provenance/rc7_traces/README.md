# rc7 micrograph traces (re-run for RC8) — the generator's masks on the papers' micrographs

opencif 1.0.0rc8, code stamp `593d977fc83ec217`; the masks are the shipped RC8 DXFs (`runs/`, sha256 in every record) and the Fig. 4 variants' builds (`variants/runs/`: the rejected readings `gifford2014_fig4_pathorder` (B) and `gifford2014_fig4_start51` (C), and the folded-head figure `gifford2014_fig4_fold20`).  Nothing in this folder is hashed code.  The method is rc6's (`provenance/rc6_traces/README.md`: what a trace is, the validators V1–V5, V11, V12, the rasters and their hashes), with rc7's and RC8's changes below.  The figures in `results/` show the papers' micrographs: they belong with the dissertation's evidence, not in a public repository (`tools/package.py` puts them in the traces zips).

## What changed in RC8

- **A7: every sample counted.**  Each record carries the outline metric's coverage (`metric.coverage`: the share of the eligible samples, a mask edge inside the micrograph and clear of the annotations, on which the metric found an image outline within its window; `metric.n_eligible`) and the censored null ratios (`metric_null.censored_ratios`: the nulls' medians over every eligible sample, one that found nothing counted at the window); every caption states both, with the chamfer null ratios.  A claim needs coverage ≥ 0.5, and the claims that use the metric nulls are judged on the censored ratios (`censored_null_ratio_min`).  (Review 1 of rc7's patch, verified; the rc6 adjudication's A7, open since rc6.)
- **D5: Xia Fig. 3d (both masks) and 3g (the depicted mask) are recorded, not claimed**: their metric covers under half of its eligible samples (d rc7 0.42, d pp 0.43, g_outlet rc7 0.48) and their censored null ratios are near 1 (1.00, 1.00, 0.94): the metric cannot tell the registration from a shifted one there.  Their registrations (V1, V5) are still reproduced and checked.
- **The Gifford Fig. 4 reproduction is the nested reading** (`gifford2014_fig4_aligned`, D1 of the adjudication of the rc7 reviews): rc7's reading B (legs in path order) is the variant `gifford2014_fig4_pathorder`, built, traced and shown, not anchored except for V13.  The Fig. 4 caption reads the neighbouring leg's walls from the record: on the nested reading leg 2's side channels are +17.5 to +18.1 µm per side wider on the image (register GIF4-1, open).
- **V13, the taper-direction sign test** (`code/taper_sign.py`, run by `run_all_rc7.sh`, anchored, and repeated by `tests/test_trace_anchors.py` layer 6): PASS when the image's outer-width slope along the leg has the mask's sign beyond 3 SE with coverage >= 0.5, >= 30 positions and a mask slope of at least 1e-4 um/um; FAIL at the opposite sign; else INCONCLUSIVE.  The SE is the larger of the Newey-West SE (Bartlett kernel, lag 100 um, beyond two row pitches) and the independent-sample SE: neighbouring positions share pixels and their residuals are autocorrelated, so the independent-sample SE alone is too small (RC8 code review, M1, N2).  A sign test: the magnitude ratio is reported, not judged.  The wall window must exceed the largest wall offset (review 3, E6; the code review, m3): +-30 um on the neighbouring leg (the nested reading's walls lie 16.3-19.6 um outside the mask's, offset_abs_range_um), +-12 um on leg 4's array part.  Required: aligned leg4 PASS, aligned neighbour PASS, fold20 leg4 PASS, pathorder leg4 PASS, pathorder neighbour FAIL, rc7 leg4 PASS, start51 leg4 PASS.

| reading | leg | window µm | verdict | positions / eligible | image slope µm/µm (z; SE the larger of Newey-West and iid) | mask slope | required |
|---|---|---|---|---|---|---|---|
| aligned | leg4 | ±12 | PASS | 80 / 80 | -0.00469 (-11.6) | -0.00800 | PASS |
| aligned | neighbour | ±30 | PASS | 130 / 130 | -0.00240 (-4.1) | -0.00257 | PASS |
| fold20 | leg4 | ±12 | PASS | 80 / 80 | -0.00471 (-11.7) | -0.00800 | PASS |
| fold20 | neighbour | ±30 | FAIL | 116 / 116 | -0.00258 (-4.2) | +0.00302 | — |
| pathorder | leg4 | ±12 | PASS | 80 / 80 | -0.00466 (-11.7) | -0.00800 | PASS |
| pathorder | neighbour | ±30 | FAIL | 117 / 117 | -0.00254 (-4.1) | +0.00301 | FAIL |
| rc7 | leg4 | ±12 | PASS | 80 / 80 | -0.00471 (-11.7) | -0.00800 | PASS |
| rc7 | neighbour | ±30 | PASS | 130 / 130 | -0.00240 (-4.1) | -0.00257 | — |
| start51 | leg4 | ±12 | PASS | 80 / 80 | -0.00472 (-11.8) | -0.00800 | PASS |
| start51 | neighbour | ±30 | PASS | 116 / 116 | -0.00258 (-4.2) | -0.00303 | — |

(`rc7` is the design-level preset `gifford2014_fig4_design`, traced under its rc7 tag.)  The image's neighbour narrows the same way in every reading (each map registers leg 4 on the same image); only the mask's sign changes, and reading B's is opposite.  The magnitude ratio is not judged: the image's slope moves by about 15 % with the sample grid.

## What changed from rc6

- **The window rule.** The nearest-outline metric's half-window is at most half the drawn gap, so a sample whose own outline is missing cannot reach the neighbouring outline across the gap (register GIF4-2: rc6's ± 8 px was 22.3 µm on Gifford's Fig. 4, wider than its 20.5 µm gap).  `panel_core.window_px(k, G)` takes the largest 0.5 px step at or below G/2 (cap 8 px); every record carries `metric_window_um` and `window_rule`, and `tests/test_trace_anchors.py` checks it.  The median and 90th percentile describe the samples that found an outline (each record's `metric.n`): a departure larger than the window drops out of them and shows grey on the residual map; the count of eligible samples was not recorded in rc7 (rc6 adjudication A7; RC8 records it, above).
- **Every shipped mask a micrograph shows**: Dinh Fig. 2B on `dinh2024_v1_depicted` (square turns, tag rc7) and `dinh2024_v1_figure_aligned` (fa); Gifford's Fig. 4 inset on `gifford2014_fig4_design` (rc7) and `gifford2014_fig4_aligned` (aligned: reading B in rc7, the nested reading in RC8); Fig. 6 on `gifford2014_fig6_design` and `_aligned` (each leg on its own map, and the four legs with their turns on one map); Xia Fig. 3 on `xia2016_rs7_depicted` (rc7) and `xia2016_rs7_aligned` (pp).  rc6's Xia variants became the shipped aligned preset; none is anchored.
- **The rejected Fig. 4 reading C** (the 51 µm start, `variants/gifford2014_fig4_start51.yaml`) is built (`code/build_variant.py`) and traced, shown, not anchored: median 0.64 / p90 1.41 µm on leg 4, its neighbour (leg 2) walls -0.08, -0.54, +0.38, -1.39 µm: the widths fit, the beads reject it (register INV-2).
- **Fig. 6's leg-7 claim end** is derived from the drawn mask (`LEG7_IMAGE_SHORT_ROWS`, `LEG7_MARGIN_ROWS`; the record's `measured.leg7_claim_end`), not fixed at 6 rows.
- **`xia_side_walls.py`** takes any record and mask on its command line and is importable (the Xia package measures with it).
- **The limits** are rc6's for every panel rc6 claimed and RC8 still claims (RC8, D5: Xia d and g on the depicted mask, rc6 claims, are recorded, their registration ranges kept and their median, p90 and null limits dropped), except: Dinh p1's scale scan across (0.975, 1.015), rc6 (0.975, 1.0): the square-turn head prefers about 1.011 (rc6's arc head 0.988); the Fig. 4 design reference's p90 3.5 µm, rc6 20 µm (rc6's 18.2 was the window, GIF4-2); Dinh p1's V2 on the chamfer nulls (its metric nulls are weak on the square-turn head); the recorded Xia f on the parallel-plate mask, its scale scan across (0.94, 0.98), rc6 (0.95, 0.98), around its 0.950.  The new masks' claims (Dinh fa, Fig. 4 aligned, Fig. 6 aligned, Xia pp c/e; Xia pp d, a claim in rc7, is recorded since RC8) were set with the results in view, tighter where the mask sits closer (Fig. 4 aligned: median ≤ 1.0, p90 ≤ 2.0 µm; Fig. 6 aligned: median ≤ 1.0, p90 ≤ 2.5 µm per leg).  They are regression limits, not a blind test; what is blind is each panel's scale, its rotation and the agreement of two independent registrations (V1).

## Results (`results/`; the limits in `trace_anchors.json`)

| panel | role | V1 along / across px | median / p90 µm | coverage | window µm | null ratio metric / censored / chamfer | V5 along / across |
|---|---|---|---|---|---|---|---|
| dinh_fig2B_p1_head__rc7 | claim | -0.87 / +0.23 | 1.40 / 2.66 | 0.916 | 8.53 | 1.17 / 1.17 / 1.70 ⚑ | 0.9799 / 1.0106 |
| dinh_fig2B_p2_transition__rc7 | claim | +0.01 / -0.14 | 0.82 / 2.55 | 0.951 | 8.53 | 2.69 / 3.66 / 2.46 | 0.9873 / 1.0057 |
| dinh_fig2B_p3_linear__rc7 | claim | +0.00 / +0.11 | 0.44 / 0.89 | 0.997 | 8.53 | 6.17 / 6.82 / 4.35 | 0.9997 / 0.9955 |
| dinh_fig2B_p4_outlet__rc7 | claim | -0.18 / -0.09 | 0.57 / 1.57 | 0.995 | 8.53 | 3.23 / 3.54 / 2.56 | 0.9958 / 1.0014 |
| dinh_fig2B_p1_head__fa | claim | -0.87 / +0.23 | 1.40 / 2.66 | 0.916 | 8.53 | 1.17 / 1.17 / 1.70 ⚑ | 0.9799 / 1.0106 |
| dinh_fig2B_p2_transition__fa | claim | +0.01 / -0.15 | 0.74 / 2.55 | 0.950 | 8.53 | 2.97 / 4.11 / 2.49 | 0.9874 / 1.0043 |
| dinh_fig2B_p3_linear__fa | claim | +0.01 / +0.10 | 0.44 / 0.89 | 0.997 | 8.53 | 6.18 / 6.79 / 4.35 | 0.9997 / 0.9955 |
| dinh_fig2B_p4_outlet__fa | claim | -0.19 / -0.10 | 0.57 / 1.58 | 0.995 | 8.53 | 3.17 / 3.46 / 2.56 | 0.9958 / 1.0017 |
| gifford_fig4_inset_leg4__rc7 | claim | -0.09 / +0.04 | 1.10 / 2.51 | 0.999 | 9.77 | 2.70 / 3.45 / 1.93 | 1.0007 / 1.0012 |
| gifford_fig4_inset_leg4__aligned | claim | -0.09 / +0.03 | 0.61 / 1.40 | 0.999 | 8.37 | 6.37 / 6.46 / 2.38 | 1.0006 / 1.0011 |
| gifford_fig6_leg7__design | claim | -0.02 / -0.36 | 1.87 / 4.26 | 0.873 | 9.05 | 1.12 / 1.36 / 1.21 ⚑ | 0.9996 / 1.0183 |
| gifford_fig6_leg5__design | claim | +0.07 / -0.01 | 1.55 / 5.47 | 0.962 | 9.05 | 1.42 / 1.57 / 1.28 ⚑ | 0.9995 / 1.0193 |
| gifford_fig6_leg3__design | claim | +0.04 / +0.01 | 1.48 / 6.50 | 0.975 | 9.05 | 1.46 / 1.55 / 1.23 ⚑ | 0.9997 / 1.0256 |
| gifford_fig6_leg1__design | claim | -0.00 / +0.10 | 1.09 / 2.33 | 0.816 | 9.05 | 2.40 / 1.99 / 1.17 ⚑ | 1.0000 / 1.0428 |
| gifford_fig6_four_legs__design | claim | -0.01 / -0.10 | 1.41 / 5.21 | 0.936 | 9.05 | 2.09 / 2.29 / 1.38 | 0.9994 / 1.0007 |
| gifford_fig6_leg7__aligned | claim | -0.00 / -0.17 | 0.64 / 1.74 | 0.872 | 6.79 | 2.22 / 4.19 / 1.63 | 0.9999 / 1.0040 |
| gifford_fig6_leg5__aligned | claim | +0.09 / -0.01 | 0.61 / 1.66 | 0.961 | 6.79 | 1.71 / 4.52 / 1.79 | 0.9995 / 1.0016 |
| gifford_fig6_leg3__aligned | claim | +0.05 / -0.03 | 0.59 / 1.64 | 0.992 | 6.79 | 2.10 / 4.28 / 1.64 | 0.9997 / 1.0022 |
| gifford_fig6_leg1__aligned | claim | +0.02 / -0.01 | 0.51 / 1.55 | 1.000 | 6.79 | 1.58 / 4.99 / 1.56 | 1.0000 / 1.0001 |
| gifford_fig6_four_legs__aligned | claim | +0.01 / -0.12 | 0.82 / 2.19 | 0.965 | 6.79 | 2.12 / 4.35 / 1.76 | 0.9991 / 0.9997 |
| xia_fig3g_outlet__rc7 | recorded | +0.55 / +1.27 | 1.73 / 3.98 | 0.476 | 7.19 | 1.04 / 0.94 / 1.05 ⚑ | 0.9980 / 1.0378 |
| xia_fig3c__rc7 | claim | +0.13 / -0.01 | 2.77 / 5.08 | 0.898 | 7.19 | 0.81 / 1.04 / 1.39 ⚑ | 0.9996 / 1.0032 |
| xia_fig3d__rc7 | recorded | +0.26 / -1.24 | 1.70 / 4.54 | 0.421 | 7.19 | 2.05 / 1.00 / 1.35 ⚑ | 1.0054 / 1.0088 |
| xia_fig3e__rc7 | recorded | +0.22 / -3.67 | 1.99 / 4.67 | 0.774 | 7.19 | 1.18 / 1.23 / 1.12 ⚑ | 1.0000 / 1.0078 |
| xia_fig3f__rc7 | recorded | -0.06 / +1.18 | 2.40 / 5.69 | 0.680 | 7.98 | 1.07 / 1.12 / 1.02 ⚑ | 1.0034 / 0.9628 |
| xia_fig3g_outlet__pp | recorded | +0.23 / +1.83 | 1.26 / 3.57 | 0.511 | 7.19 | 0.91 / 1.15 / 1.32 ⚑ | 0.9969 / 1.0043 |
| xia_fig3c__pp | claim | +0.13 / -0.01 | 2.78 / 4.99 | 0.903 | 7.19 | 0.81 / 1.04 / 1.38 ⚑ | 0.9996 / 1.0034 |
| xia_fig3d__pp | recorded | -0.25 / +0.39 | 2.99 / 4.76 | 0.435 | 7.19 | 1.06 / 1.00 / 1.37 ⚑ | 1.0027 / 0.9988 |
| xia_fig3e__pp | claim | +0.19 / +0.94 | 2.02 / 4.11 | 0.897 | 7.19 | 0.53 / 1.18 / 1.21 ⚑ | 1.0000 / 0.9981 |
| xia_fig3f__pp | recorded | -0.06 / +1.17 | 2.41 / 5.71 | 0.681 | 7.98 | 1.07 / 1.13 / 1.02 ⚑ | 1.0032 / 0.9500 |

⚑ a null ratio below 1.3 (package check, finding 7).  `null_flags` in `trace_anchors.json` lists each, with whether the panel's V2 claim uses that null: dinh_fig2B_p1_head__rc7 metric 1.17; dinh_fig2B_p1_head__rc7 censored 1.17; dinh_fig2B_p1_head__fa metric 1.17; dinh_fig2B_p1_head__fa censored 1.17; gifford_fig6_leg7__design chamfer 1.21 (V2); gifford_fig6_leg7__design metric 1.12; gifford_fig6_leg5__design chamfer 1.28 (V2); gifford_fig6_leg3__design chamfer 1.23 (V2); gifford_fig6_leg1__design chamfer 1.17 (V2); xia_fig3g_outlet__rc7 chamfer 1.05; xia_fig3g_outlet__rc7 metric 1.04; xia_fig3g_outlet__rc7 censored 0.94; xia_fig3c__rc7 metric 0.81; xia_fig3c__rc7 censored 1.04; xia_fig3d__rc7 censored 1.00; xia_fig3e__rc7 chamfer 1.12; xia_fig3e__rc7 metric 1.18; xia_fig3e__rc7 censored 1.23; xia_fig3f__rc7 chamfer 1.02; xia_fig3f__rc7 metric 1.07; xia_fig3f__rc7 censored 1.12; xia_fig3g_outlet__pp metric 0.91; xia_fig3g_outlet__pp censored 1.15; xia_fig3c__pp metric 0.81; xia_fig3c__pp censored 1.04; xia_fig3d__pp metric 1.06; xia_fig3d__pp censored 1.00; xia_fig3e__pp chamfer 1.21 (V2); xia_fig3e__pp metric 0.53; xia_fig3e__pp censored 1.18; xia_fig3f__pp chamfer 1.02; xia_fig3f__pp metric 1.07; xia_fig3f__pp censored 1.13.

**Xia Fig. 3a** (V12, each traced panel against its box): rc7: c, d, e, f inside their boxes; outside: g 8.20 of its box's half-size; pp: c, d, e, f inside their boxes; outside: g 1.09 of its box's half-size.

**Xia panels f and g, the side channels' outer walls** (`results/xia/xia_side_walls.json`, image inside the mask's, on each panel's own registration; register INV-4, open): f_rc7 left +23.6 µm, right +19.7 µm (1.774 µm/px); g_rc7 left -26.4 µm (1.798 µm/px); f_pp left +39.9 µm, right +36.0 µm (1.774 µm/px); g_pp left -0.3 µm (1.798 µm/px).

## Reproducing

```
export OPENCIF_PAPER_RASTERS=<the folder of dinh5-000.png, gif5-001.png, gif7-000.png, xia5-000.png>   # rc6 README
python3 tests/test_trace_anchors.py          # every registration repeated from the rasters and the DXFs; V13 (layer 6)
```

To regenerate the records: `provenance/rc7_traces/code/run_all_rc7.sh RASTERS` (resumable by its logs; it builds the Fig. 4 variants first and runs V13 last), then `RASTERS=… python3 provenance/rc7_traces/code/make_anchors.py . provenance/rc7_traces/results provenance/rc7_traces/trace_anchors.json` and `python3 provenance/rc7_traces/code/readme_rc7.py .` (this file).

Licences: as rc6 (`provenance/rc6_traces/README.md`): Xia 2016 CC BY 4.0; Dinh 2024 and Gifford 2014 © The Royal Society of Chemistry.  The rasters are not redistributed here.
