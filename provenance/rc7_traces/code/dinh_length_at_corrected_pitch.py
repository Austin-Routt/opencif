"""The Dinh length paragraph's numbers at the row pitch measured with the bars read at their anti-aliased length
(45.1 +- 0.2 um on Fig. 2B's panels 2, 3 and 4: 45.16, 45.09, 45.08) instead of rc5's 46.0 +- 1.0 um.  The premises
and the functions are rc5's (provenance/rc5_scripts/dinh_taper_alternatives.py, dinh_length_bound.py): the frozen
cutoff law's cubic kernel at the plateau prefactor, the lamina the same at every gap, fixed by the printed
f* = 1.04e-4 at 75 um; every centre width >= 52 um; two gaps per row; every gap row at least the pitch long.  Rows do
not depend on the pitch; lengths scale with it, and the share of rows that may taper depends on it through the number
of rows that fit in 75.8 mm.  Also Reviewer A's head corollary (E2) at the corrected pitch.

usage: python3 dinh_length_at_corrected_pitch.py OPENCIF_ROOT OUT.json"""
import json
import math
import os
import sys

R, OUT = sys.argv[1:3]
sys.path[:0] = [R]
from opencif.physics import schedule as S

fstar, w_in, w_min = 1.04e-4, 75.0, 52.0
x_star = S.x_of(fstar)


def f_at(w):
    x = x_star * w_in / w
    return 3 * x ** 2 - 2 * x ** 3


f_max = f_at(w_min)


def F_of(widths, extra=0.0):
    return 1.0 - math.exp(sum(math.log(1.0 - 2.0 * f_at(w)) for w in widths) + extra * math.log(1.0 - 2.0 * f_max))


def linear_all(N):
    return [w_in - (w_in - w_min) * i / (N - 1) for i in range(N)]


def linear_then_flat(N, q):
    n = max(2, int(round(q * N)))
    return [w_in - (w_in - w_min) * i / (n - 1) if i < n else w_min for i in range(N)]


def rows_needed(F, shape):
    lo, hi = 2, 20000
    while F_of(shape(hi)) < F:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if F_of(shape(mid)) >= F else (mid, hi)
    return hi


def q_fit(F, n_rows, credit=0.0):
    for k in range(1000, -1, -1):
        if F_of(linear_then_flat(n_rows, k / 1000.0), credit) >= F:
            return k / 1000.0
    return None


credits = {}
for k in ('depicted', 'specified'):
    jp = os.path.join(R, 'provenance', 'rc5_scripts', 'dinh_length_bound_%s.json' % k)
    if os.path.exists(jp):
        credits[k] = float(json.load(open(jp))['rows_credited_to_this_head'])
printed = 75.8
res = dict(f_max=f_max, printed_mm=printed, head_credit_rows=credits, pitches={})
for pt in (46.0, 45.3, 45.1, 44.9):
    rows = []
    for lab, F in (('presets (outlet identity)', 0.5592819159228324), ('constant-shear taper reaching 52 um', 0.5105),
                   ('measured mean', 0.505), ('measured mean minus one SD', 0.472)):
        n_b = rows_needed(F, lambda N: [w_min] * N)
        n_l = rows_needed(F, linear_all)
        n_fit = int(printed * 1000 / pt)
        row = dict(F_reading=lab, F=round(F, 4), rows_every_row=n_b, bound_mm=round(n_b * pt / 1000, 2),
                   bound_mm_credited={k: round((n_b - c) * pt / 1000, 2) for k, c in credits.items()},
                   rows_linear=n_l, linear_taper_mm=round(n_l * pt / 1000, 2), rows_in_printed_length=n_fit,
                   largest_share_tapering=q_fit(F, n_fit),
                   largest_share_tapering_with_head_credit={k: q_fit(F, n_fit, c) for k, c in credits.items()})
        rows.append(row)
    # Reviewer A's corollary (E2): the reconstructed head counted at its drawn length, the rest of 75.8 mm at f_max
    head = dict(depicted=dict(length_mm=25.676, P_after_head=0.944), specified=dict(length_mm=28.614, P_after_head=0.944))
    cor = {}
    for k, h in head.items():
        n_rest = int((printed - h['length_mm']) * 1000 / pt)
        cor[k] = dict(rows_after_head=n_rest, F_max=round(1.0 - h['P_after_head'] * (1.0 - 2.0 * f_max) ** n_rest, 4))
    res['pitches']['%.1f' % pt] = dict(readings=rows, head_corollary_E2=cor)
json.dump(res, open(OUT, 'w'), indent=1)
for pt, d in res['pitches'].items():
    print('pitch %s um' % pt)
    for r in d['readings']:
        print('  F %.4f  bound %.2f mm (credited %s)  linear %.2f mm  share tapering %s (credited %s)' % (
            r['F'], r['bound_mm'], r['bound_mm_credited'], r['linear_taper_mm'], r['largest_share_tapering'],
            r['largest_share_tapering_with_head_credit']))
    print('  E2 head corollary', d['head_corollary_E2'])
