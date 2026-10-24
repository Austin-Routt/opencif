"""The model band, rc4 vs rc5, per drawn design: total rows, rows among the last 100, largest change and its row,
rows losing their cutoff, share of F; with the outlet's stretch-end term.  The rows come from the records' lists
(validation.model_scenarios junction_band_lo/hi); records before code 56e1980d8a064b35 list at most 40 per scenario,
so on those the count of last rows can fall short of the band's count.
Usage: python3 band_rc5.py RC4_RUNS RC5_RUNS"""
import json, glob, os, sys
def recs(root):
    out = {}
    for p in sorted(glob.glob(os.path.join(root, '*', '*.json'))):
        n = os.path.basename(os.path.dirname(p))
        if os.path.basename(p) != n + '.json': continue
        r = json.load(open(p))
        if 'validation' not in r or 'element' in r: continue
        out[n] = r
    return out
A = recs(sys.argv[1] if len(sys.argv) > 2 else 'opencif_rc4/runs')
B = recs(sys.argv[2] if len(sys.argv) > 2 else 'opencif_rc5/runs')
def band(r):
    ms = r['validation']['model_scenarios']; N = r['design']['N']
    lo, hi = ms['junction_band_lo'], ms['junction_band_hi']
    rows = sorted(set(lo['newly_outside_rows']) | set(hi['newly_outside_rows']))
    n = max(lo['rows_newly_outside_tol'], hi['rows_newly_outside_tol'])
    mx = max(lo['max_cutoff_change'] or 0, hi['max_cutoff_change'] or 0)
    mrow = lo['max_change_row'] if (lo['max_cutoff_change'] or 0) >= (hi['max_cutoff_change'] or 0) else hi['max_change_row']
    lost = max(lo.get('rows_losing_their_cutoff', 0), hi.get('rows_losing_their_cutoff', 0))
    sh = max(lo.get('share_of_F_in_those_rows') or 0, hi.get('share_of_F_in_those_rows') or 0)
    last = [x for x in rows if x > N - 100]
    return n, len(last), mx, mrow, N, lost, sh, rows
print('%-46s %-17s %5s | rc4 n/max | rc5 n (last) max@row lost shF' % ('design', 'outlet end', 'c'))
summ = []
for n in sorted(B):
    se = (B[n]['design'].get('outlet') or {}).get('stretch_end') or {}
    a, b = band(A[n]), band(B[n])
    seam_same = [x for x in a[7]] == [x for x in b[7] if x <= b[4] - 100]
    summ.append((n, se.get('kind', 'none'), se.get('c', 0.0), a, b, seam_same))
    print('%-46s %-17s %5.3f | %3d %5.3f | %3d (%3d) %5.3f@%-5d %2d %.4f %s' % (n, se.get('kind', 'none'), se.get('c', 0.0), a[0], a[2], b[0], b[1], b[2], b[3], b[5], b[6], '' if seam_same else 'SEAM-CHANGED'))
step = [s for s in summ if s[1] != 'none']
nost = [s for s in summ if s[1] == 'none']
print('designs with a stretch-end term: %d; band names last rows on %d of them' % (len(step), sum(1 for s in step if s[4][1])))
print('  last rows named: %d-%d; band max over those designs (any row) %.3f-%.3f; lost up to %d' % (
    min(s[4][1] for s in step), max(s[4][1] for s in step), min(s[4][2] for s in step), max(s[4][2] for s in step), max(s[4][5] for s in step)))
print('designs without: %d; band names last rows on %d (lead arcs)' % (len(nost), sum(1 for s in nost if s[4][1])))
lr = [s for s in nost if s[4][1]]
if lr:
    print('  ', [(s[0], s[4][1], round(s[4][2], 4), s[4][3]) for s in lr])
print('seam rows unchanged on all:', all(s[5] for s in summ), [s[0] for s in summ if not s[5]])
