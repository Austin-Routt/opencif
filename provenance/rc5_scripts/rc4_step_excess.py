"""rc4's unmodelled outlet-stretch end, design by design: the width change w_s(N) -> lead width read from each rc4 run
pickle, its excess from the rc5 coefficient (cifsched.outlet_exp_c, imported from the rc5 tree), as a share of the
side outlet's resistance, with the rc4 record's own d(D_N/D)/d ln R_s,out (a linearised size: beyond a few % the last
rows lose their cutoff).  Usage: python3 rc4_step_excess.py <rc4 tree> <rc5 tree>"""
import sys, os, glob, json, pickle
RC4, RC5 = sys.argv[1], sys.argv[2]
sys.path[:0] = [RC5, os.path.join(RC5, 'opencif'), os.path.join(RC5, 'opencif', 'engine')]
from opencif.engine import cifsched as S
rows = []
for pk in sorted(glob.glob(os.path.join(RC4, 'runs', '*', '*.pkl'))):
    n = os.path.basename(os.path.dirname(pk))
    rec = os.path.join(RC4, 'runs', n, n + '.json')
    if not os.path.exists(rec):
        continue
    r = json.load(open(rec))
    if 'validation' not in r:
        continue
    a = pickle.load(open(pk, 'rb'))['arrays']
    o, spec, N = a['outlet'], a['spec'], int(a['N'])
    style = getattr(spec, 'outlet_style', 'teardrop') if o.get('kind') == 'combined' else 'separate'
    w_last, w_lead = float(a['ws'][N - 2]), float(o.get('w_lead') or 0.0)
    c, deep, shal, lam, kind = S.outlet_exp_c(spec, w_last, w_lead, parts=True)
    wn = min(w_last, w_lead)
    dR = c * wn * S.rr(wn, spec.h) * S.UM
    Rs = r['validation']['outlet_measured']['RsA_out']
    der = (r['validation'].get('model_scenarios') or {}).get('dDlast_dlnRside_out')
    rows.append(dict(name=n, outlet=style, w_side=round(w_last, 2), w_lead=round(w_lead, 1), ratio=round(max(w_last, w_lead) / wn, 3),
                     kind=kind, c=round(float(c), 4), c_deep=round(float(deep), 4), c_shallow=round(float(shal), 4),
                     depth_weight=round(float(lam), 3), share_of_Rs_out_pct=round(100 * dR / Rs, 3),
                     dDN_dlnR=round(der, 1) if der is not None else None,
                     linearised_last_row_change_D=round(der * dR / Rs, 2) if der is not None else None))
json.dump(rows, open('rc4_step_excess.json', 'w'), indent=1)
print('| design | outlet | w_s(N) | lead | ratio | change | c (narrow w) | deep / shallow / weight | % of R_s,out | dD_N/d ln R | linearised (D) |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for x in rows:
    if x['c'] > 0:
        print('| %s | %s | %.2f | %.1f | %.2f | %s | %.3f | %.3f / %.3f / %.3f | %.2f | %s | %s |' % (
            x['name'], x['outlet'], x['w_side'], x['w_lead'], x['ratio'], x['kind'], x['c'], x['c_deep'], x['c_shallow'],
            x['depth_weight'], x['share_of_Rs_out_pct'], x['dDN_dlnR'], x['linearised_last_row_change_D']))
st = [x for x in rows if x['c'] > 0]
print('\n%d of %d drawn designs change width at the stretch end (every teardrop and ring outlet; no separate outlet, all '
      'with w_s(N) >= 40 um); excess %.2f-%.2f %% of R_s,out' % (len(st), len(rows), min(x['share_of_Rs_out_pct'] for x in st),
                                                              max(x['share_of_Rs_out_pct'] for x in st)))
