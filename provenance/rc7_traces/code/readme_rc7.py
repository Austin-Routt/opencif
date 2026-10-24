"""Write provenance/rc7_traces/README.md from trace_anchors.json and the records in results/ (rc7; RC8).

usage: python3 -I readme_rc7.py ROOT"""
import json
import os
import sys

ROOT = sys.argv[1]
TR = os.path.join(ROOT, 'provenance', 'rc7_traces')
A = json.load(open(os.path.join(TR, 'trace_anchors.json')))
RES = os.path.join(TR, 'results')
sw = json.load(open(os.path.join(RES, 'xia', 'xia_side_walls.json')))['panels']
s51 = json.load(open(os.path.join(RES, 'gifford', 'gifford_results_start51.json')))
po = json.load(open(os.path.join(RES, 'gifford', 'gifford_results_pathorder.json')))


def recs():
    out = {}
    for d, fs in (('dinh', ['dinh_results_rc7.json', 'dinh_results_fa.json']),
                  ('gifford', ['gifford_results_rc7.json', 'gifford_results_aligned.json', 'gifford_results_start51.json',
                               'gifford_results_pathorder.json', 'gifford_results_fold20.json',
                               'gifford6_results_design.json', 'gifford6_results_design_onemap.json',
                               'gifford6_results_aligned.json', 'gifford6_results_aligned_onemap.json']),
                  ('xia', ['xia_results_%s.json' % t for t in ('rc7', 'pp')] +
                   ['xia_results_fig3%s_%s.json' % (L, t) for t in ('rc7', 'pp') for L in 'cdef'])):
        for f in fs:
            r = json.load(open(os.path.join(RES, d, f)))
            for v in ([r] if 'name' in r else list(r.values())):
                if isinstance(v, dict) and 'name' in v:
                    out[v['name']] = v
    return out


R = recs()
L = []
w = L.append
w('# rc7 micrograph traces (re-run for RC8) — the generator\'s masks on the papers\' micrographs')
w('')
w('opencif 1.0.0rc8, code stamp `%s`; the masks are the shipped RC8 DXFs (`runs/`, sha256 in every record) and the '
  'Fig. 4 variants\' builds (`variants/runs/`: the rejected readings `gifford2014_fig4_pathorder` (B) and '
  '`gifford2014_fig4_start51` (C), and the folded-head figure `gifford2014_fig4_fold20`).  Nothing in this folder is '
  'hashed code.  The method is rc6\'s (`provenance/rc6_traces/README.md`: what a trace is, the validators V1–V5, V11, '
  'V12, the rasters and their hashes), with rc7\'s and RC8\'s changes below.  The figures in `results/` show the papers\' '
  'micrographs: they belong with the dissertation\'s evidence, not in a public repository (`tools/package.py` puts them '
  'in the traces zips).' % A['code_stamp']['sha256'])
w('')
w('## What changed in RC8')
w('')
w('- **A7: every sample counted.**  Each record carries the outline metric\'s coverage (`metric.coverage`: the share of '
  'the eligible samples, a mask edge inside the micrograph and clear of the annotations, on which the metric found an '
  'image outline within its window; `metric.n_eligible`) and the censored null ratios (`metric_null.censored_ratios`: '
  'the nulls\' medians over every eligible sample, one that found nothing counted at the window); every caption states '
  'both, with the chamfer null ratios.  A claim needs coverage ≥ %.1f, and the claims that use the metric nulls are '
  'judged on the censored ratios (`censored_null_ratio_min`).  (Review 1 of rc7\'s patch, verified; the rc6 '
  'adjudication\'s A7, open since rc6.)' % 0.5)
w('- **D5: Xia Fig. 3d (both masks) and 3g (the depicted mask) are recorded, not claimed**: their metric covers under '
  'half of its eligible samples (%s) and their censored null ratios are near 1 (%s): the metric cannot tell the '
  'registration from a shifted one there.  Their registrations (V1, V5) are still reproduced and checked.' % (
      ', '.join('%s %.2f' % (n.replace('xia_fig3', '').replace('__', ' '), A['panels'][n]['anchored']['coverage'])
                for n in ('xia_fig3d__rc7', 'xia_fig3d__pp', 'xia_fig3g_outlet__rc7')),
      ', '.join('%.2f' % A['panels'][n]['null_ratios'].get('censored', float('nan'))
                for n in ('xia_fig3d__rc7', 'xia_fig3d__pp', 'xia_fig3g_outlet__rc7'))))
w('- **The Gifford Fig. 4 reproduction is the nested reading** (`gifford2014_fig4_aligned`, D1 of the adjudication of '
  'the rc7 reviews): rc7\'s reading B (legs in path order) is the variant `gifford2014_fig4_pathorder`, built, traced and '
  'shown, not anchored except for V13.  The Fig. 4 caption reads the neighbouring leg\'s walls from the record: on the '
  'nested reading leg 2\'s side channels are %s µm per side wider on the image (register GIF4-1, open).' % (
      ' to '.join('%+.1f' % v for v in (lambda b: (b[0], b[-1]))(sorted(
          x['mean_signed_um'] for x in R['gifford_fig4_inset_leg4__aligned']['measured']['leg2_side_walls']['walls'].values()
          if x['n'] >= 50)))))
w('- **V13, the taper-direction sign test** (`code/taper_sign.py`, run by `run_all_rc7.sh`, anchored, and repeated by '
  '`tests/test_trace_anchors.py` layer 6): ' + A['v13_note'].split(': ', 1)[1])
w('')
w('| reading | leg | window µm | verdict | positions / eligible | image slope µm/µm (z; SE the larger of Newey-West and iid) | mask slope | required |')
w('|---|---|---|---|---|---|---|---|')
for k_, v_ in sorted(A['v13'].items()):
    w('| %s | %s | ±%g | %s | %d / %d | %s | %s | %s |' % (
        k_.rsplit('_', 1)[0], v_['leg'], v_['window_um'], v_['verdict'], v_['positions'], v_['eligible'],
        ('%+.5f (%.1f)' % (v_['slope_image'], v_['z_image'])) if v_.get('slope_image') is not None else '—',
        ('%+.5f' % v_['slope_mask']) if v_.get('slope_mask') is not None else '—', v_.get('required') or '—'))
w('')
w('(`rc7` is the design-level preset `gifford2014_fig4_design`, traced under its rc7 tag.)  The image\'s neighbour narrows '
  'the same way in every reading (each map registers leg 4 on the same image); only the mask\'s sign changes, and reading B\'s is opposite.  The magnitude ratio is not '
  'judged: the image\'s slope moves by about 15 % with the sample grid.')
w('')
w('## What changed from rc6')
w('')
w('- **The window rule.** The nearest-outline metric\'s half-window is at most half the drawn gap, so a sample whose own '
  'outline is missing cannot reach the neighbouring outline across the gap (register GIF4-2: rc6\'s ± 8 px was 22.3 µm on '
  'Gifford\'s Fig. 4, wider than its 20.5 µm gap).  `panel_core.window_px(k, G)` takes the largest 0.5 px step at or below '
  'G/2 (cap 8 px); every record carries `metric_window_um` and `window_rule`, and `tests/test_trace_anchors.py` checks it.  The median '
  'and 90th percentile describe the samples that found an outline (each record\'s `metric.n`): a departure larger than '
  'the window drops out of them and shows grey on the residual map; the count of eligible samples was not recorded in '
  'rc7 (rc6 adjudication A7; RC8 records it, above).')
w('- **Every shipped mask a micrograph shows**: Dinh Fig. 2B on `dinh2024_v1_depicted` (square turns, tag rc7) and '
  '`dinh2024_v1_figure_aligned` (fa); Gifford\'s Fig. 4 inset on `gifford2014_fig4_design` (rc7) and `gifford2014_fig4_aligned` '
  '(aligned: reading B in rc7, the nested reading in RC8); Fig. 6 on `gifford2014_fig6_design` and `_aligned` (each leg on its own map, and the four legs '
  'with their turns on one map); Xia Fig. 3 on `xia2016_rs7_depicted` (rc7) and `xia2016_rs7_aligned` (pp).  rc6\'s Xia '
  'variants became the shipped aligned preset; none is anchored.')
w('- **The rejected Fig. 4 reading C** (the 51 µm start, `variants/gifford2014_fig4_start51.yaml`) is built '
  '(`code/build_variant.py`) and traced, shown, not anchored: median %.2f / p90 %.2f µm on leg 4, its neighbour (leg 2) '
  'walls %s µm: the widths fit, the beads reject it (register INV-2).' % (
      s51['metric']['median_abs_um'], s51['metric']['p90_abs_um'],
      ', '.join('%+.2f' % v['mean_signed_um'] for v in s51['measured']['leg2_side_walls']['walls'].values())))
w('- **Fig. 6\'s leg-7 claim end** is derived from the drawn mask (`LEG7_IMAGE_SHORT_ROWS`, `LEG7_MARGIN_ROWS`; the record\'s '
  '`measured.leg7_claim_end`), not fixed at 6 rows.')
w('- **`xia_side_walls.py`** takes any record and mask on its command line and is importable (the Xia package measures '
  'with it).')
w('- **The limits** are rc6\'s for every panel rc6 claimed and RC8 still claims (RC8, D5: Xia d and g on the depicted '
  'mask, rc6 claims, are recorded, their registration ranges kept and their median, p90 and null limits dropped), except: Dinh p1\'s scale scan across (0.975, 1.015), rc6 (0.975, '
  '1.0): the square-turn head prefers about 1.011 (rc6\'s arc head 0.988); the Fig. 4 design reference\'s p90 3.5 µm, rc6 '
  '20 µm (rc6\'s 18.2 was the window, GIF4-2); Dinh p1\'s V2 on the chamfer nulls (its metric nulls are weak on the '
  'square-turn head); the recorded Xia f on the parallel-plate mask, its scale scan across (0.94, 0.98), rc6 (0.95, '
  '0.98), around its 0.950.  The new masks\' claims (Dinh fa, Fig. 4 aligned, Fig. 6 aligned, Xia pp c/e; Xia pp d, a claim in rc7, is recorded since RC8) were set with the '
  'results in view, tighter where the mask sits closer (Fig. 4 aligned: median ≤ 1.0, p90 ≤ 2.0 µm; Fig. 6 aligned: '
  'median ≤ 1.0, p90 ≤ 2.5 µm per leg).  They are regression limits, not a blind test; what is blind is each panel\'s '
  'scale, its rotation and the agreement of two independent registrations (V1).')
w('')
w('## Results (`results/`; the limits in `trace_anchors.json`)')
w('')
w('| panel | role | V1 along / across px | median / p90 µm | coverage | window µm | null ratio metric / censored / chamfer | V5 along / across |')
w('|---|---|---|---|---|---|---|---|')
for n, P in A['panels'].items():
    r = R.get(n)
    a, nr = P['anchored'], P['null_ratios']
    v5 = a.get('V5') or [None, None]
    flag = '' if min(nr.values()) >= 1.3 else ' ⚑'
    w('| %s | %s | %+.2f / %+.2f | %.2f / %.2f | %s | %.2f | %s / %s / %.2f%s | %s |' % (
        n, P['role'], a['V1_along_px'], a['V1_across_px'], a['median_abs_um'], a['p90_abs_um'],
        ('%.3f' % a['coverage']) if a.get('coverage') is not None else '—', a['metric_window_um'] or 0.0,
        ('%.2f' % nr['metric']) if 'metric' in nr else '—', ('%.2f' % nr['censored']) if 'censored' in nr else '—',
        nr['chamfer'], flag, ('%.4f / %.4f' % tuple(v5)) if v5[0] is not None else '—'))
w('')
w('⚑ a null ratio below 1.3 (package check, finding 7).  `null_flags` in `trace_anchors.json` lists each, with whether the '
  'panel\'s V2 claim uses that null: %s.' % ('; '.join('%s %s %.2f%s' % (f['panel'], f['null'], f['ratio'], ' (V2)' if f['used_for_V2'] else '')
                                                    for f in A['null_flags'])))
w('')
w('**Xia Fig. 3a** (V12, each traced panel against its box): ' + '; '.join(
    '%s: %s' % (t, v['V12_expected']) for t, v in A['fig3a'].items()) + '.')
w('')
w('**Xia panels f and g, the side channels\' outer walls** (`results/xia/xia_side_walls.json`, image inside the mask\'s, on '
  'each panel\'s own registration; register INV-4, open): ' + '; '.join(
      '%s %s (%.3f µm/px)' % (k, ', '.join('%s %+.1f µm' % (x['side'], x['image_inside_um']) for x in v['walls']), v['um_per_px'])
      for k, v in sw.items()) + '.')
w('')
w('## Reproducing')
w('')
w('```')
w('export OPENCIF_PAPER_RASTERS=<the folder of dinh5-000.png, gif5-001.png, gif7-000.png, xia5-000.png>   # rc6 README')
w('python3 tests/test_trace_anchors.py          # every registration repeated from the rasters and the DXFs; V13 (layer 6)')
w('```')
w('')
w('To regenerate the records: `provenance/rc7_traces/code/run_all_rc7.sh RASTERS` (resumable by its logs; it builds the '
  'Fig. 4 variants first and runs V13 last), then `RASTERS=… python3 provenance/rc7_traces/code/make_anchors.py . provenance/rc7_traces/results '
  'provenance/rc7_traces/trace_anchors.json` and `python3 provenance/rc7_traces/code/readme_rc7.py .` (this file).')
w('')
w('Licences: as rc6 (`provenance/rc6_traces/README.md`): Xia 2016 CC BY 4.0; Dinh 2024 and Gifford 2014 © The Royal '
  'Society of Chemistry.  The rasters are not redistributed here.')
open(os.path.join(TR, 'README.md'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
print('wrote', os.path.join(TR, 'README.md'))
