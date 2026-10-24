"""rc5 vs rc6 record statistics quoted in the docs (outlet balance, band split, AR, robustness, lengths).
Usage: python3 stats_rc6.py RC5_RUNS RC6_RUNS
(the band's rows come from validation.model_scenarios junction_band_lo/hi; records before code 56e1980d8a064b35 list
at most 40 rows per scenario)"""
import json, glob, os, sys
def recs(root):
    out = {}
    for p in sorted(glob.glob(os.path.join(root, '*', '*.json'))):
        n = os.path.basename(os.path.dirname(p))
        if os.path.basename(p) != n + '.json':
            continue
        r = json.load(open(p))
        if 'element' in r or 'refused' in r or 'validation' not in r:
            continue
        out[n] = r
    return out
def chk(r, pre):
    for c in r['assessment']['checks']:
        if c['check'].startswith(pre):
            return c
RA = sys.argv[1]
RB = sys.argv[2]
for lab, root in (('rc5', RA), ('rc6', RB)):
    R = recs(root)
    print('=====', lab, len(R), 'drawn records')
    der, rows, dmax, fsh, lastw, lost = [], [], [], [], [], []
    seam_rows, last_rows, seam_max, last_max, sh_band = {}, {}, {}, {}, {}
    for n, r in R.items():
        V = r['validation']; d = r['design']; N, NM = d['N'], d['NM']
        ms = V['model_scenarios']
        der.append((ms['dDlast_dlnRside_out'], n))
        p, m = ms['outlet_side_resistance_+1%'], ms['outlet_side_resistance_-1%']
        rows.append((max(p['rows_newly_outside_tol'], m['rows_newly_outside_tol']), n))
        dmax.append((max(p['max_cutoff_change'] or 0, m['max_cutoff_change'] or 0), n))
        fsh.append((max(p['share_of_F_in_those_rows'] or 0, m['share_of_F_in_those_rows'] or 0), n))
        allr = sorted(set(p['newly_outside_rows']) | set(m['newly_outside_rows']))
        if allr: lastw.append((N - min(allr) + 1, n))
        lost.append((max(p.get('rows_losing_their_cutoff', 0), m.get('rows_losing_their_cutoff', 0)), n))
        b = chk(r, '1-D model band')['value']
        jb = [v for k, v in b.items() if k.startswith('junction band')]
        if jb:
            jb = jb[0]
            rr = sorted(set(ms['junction_band_lo']['newly_outside_rows']) | set(ms['junction_band_hi']['newly_outside_rows']))
            lastr = [x for x in rr if x > N - 100]
            seamr = [x for x in rr if x <= N - 100]
            last_rows[n] = lastr; seam_rows[n] = seamr
            last_max[n] = jb['max_cutoff_change_over_D']; sh_band[n] = jb.get('share_of_F_in_those_rows')
    rng = lambda L: (min(L)[0], min(L)[1], max(L)[0], max(L)[1])
    print('dDlast_dlnRside_out range', rng(der))
    print('outlet +-1% rows newly outside', rng(rows))
    print('outlet +-1% max change D', rng(dmax))
    nx = [x for x in dmax if not x[1].startswith('xia')]
    print('   without Xia', rng(nx), ' Xia', [x for x in dmax if x[1].startswith('xia')])
    print('outlet +-1% share of F', rng(fsh))
    print('within the last', rng(lastw), ' non-Xia', rng([x for x in lastw if not x[1].startswith('xia')]))
    print('rows losing cutoff', rng(lost))
    # band split
    withlast = {n: v for n, v in last_rows.items() if v}
    print('band names last rows on %d designs' % len(withlast))
    if withlast:
        print('   last-row counts', sorted((len(v), n) for n, v in withlast.items())[:3], '...', sorted((len(v), n) for n, v in withlast.items())[-3:])
    # AR check
    ar = {n: chk(r, 'aspect ratio of the narrowest')['value'] for n, r in R.items()}
    print('AR at:', {k: sum(1 for v in ar.values() if v['at'] == k) for k in sorted(set(v['at'] for v in ar.values()))})
    over10 = sorted((v['depth_over_width'], n) for n, v in ar.items() if v['depth_over_width'] > 10)
    print('AR > 10:', over10)
    # robustness ranges
    rb = {n: r['validation']['robustness_design_network'] for n, r in R.items()}
    # lengths of Dinh devices
    for n in ('dinh2024_v1_depicted', 'dinh2024_v1_specified', 'dinh_targets_sealed_head'):
        if n in R:
            print(n, 'footprint', [round(x, 3) for x in R[n]['validation']['drawn_parameters']['footprint_mm']])
    # Re bound / dean
    g = [(r['validation']['gates'].get('Re_any_channel_bound'), n) for n, r in R.items()]
    print('Re_any <= 10 on', sum(1 for x in g if x[0] is not None and x[0] <= 10), 'of', len(g))
