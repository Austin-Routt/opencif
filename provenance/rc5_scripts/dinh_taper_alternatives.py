"""rc5 (after the rc5 reviews): how the Dinh V1 row count depends on the centre channel's taper, with every premise of
dinh_length_bound.py kept (the frozen cutoff law's cubic kernel with the lamina the same at every gap, at the plateau
prefactor kappa_g, fixed by the printed f* = 1.04e-4 at 75 um, so kappa cancels; the march's recursion P(i+1) = P(i) (1 - 2 f(i)); every row at least the measured 46.0 um pitch) and the
caption's "narrowed from 75 um near the inlet to 52 um near the outlet" read as three taper shapes:
  (1) the bound itself: every row at 52 um (f_max), the least number of rows any taper can have;
  (2) a centre width linear in the row index from 75 um at row 1 to 52 um at the last row;
  (3) linear from 75 um to 52 um over the first q of the rows, then 52 um: the largest q whose rows still fit in the
      printed 75.8 mm if the rows had the whole length (no port, filter, head turn or outlet);
and, for any shape, the constant-equivalent per-row extraction the rows need to fit in that length (the condition is
on sum ln(1 - 2 f_i), so this is the f a constant channel would need), as a fraction of f_max and as the constant
centre width that extracts it (a width below 52 um, against the premise, means no taper fits).  Because f(w) is
convex, a taper whose MEAN width exceeds that width can still fit: the condition is on extraction, not on mean width.
No head is modelled (every row is a gap row at the pitch, at the constant lamina, so no row extracts more than f_max).
After an independent check of the rc5-review report: (3) is repeated at the pitch's +-1 um (45 and 47 um) and with
the rows the presets' reconstructed heads draw above f_max credited as extra f_max rows (dinh_length_bound_*.json,
rows_credited_to_this_head; a departure from the constant lamina).
Usage: python3 dinh_taper_alternatives.py OPENCIF_ROOT [OUT.json]"""
import json, math, os, sys
R = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else 'dinh_taper_alternatives.json'
sys.path[:0] = [R, os.path.join(R, 'opencif', 'engine')]
import cifsched as S

fstar, w_in, w_min = 1.04e-4, 75.0, 52.0
x_star = S.x_of(fstar)


def f_at(w):
    x = x_star * w_in / w
    return 3 * x ** 2 - 2 * x ** 3


f_max = f_at(w_min)


def F_of(widths, extra_fmax_rows=0.0):
    lnP = sum(math.log(1.0 - 2.0 * f_at(w)) for w in widths) + extra_fmax_rows * math.log(1.0 - 2.0 * f_max)
    return 1.0 - math.exp(lnP)


def linear_all(N):
    return [w_in - (w_in - w_min) * i / (N - 1) for i in range(N)]


def linear_then_flat(N, q):
    n = max(2, int(round(q * N)))
    return [w_in - (w_in - w_min) * i / (n - 1) if i < n else w_min for i in range(N)]


def rows_needed(F, shape):
    lo, hi = 2, 20000
    while F_of(shape(hi)) < F:
        hi *= 2
    while hi - lo > 1:                       # F_of(shape(N)) rises with N for these shapes
        mid = (lo + hi) // 2
        if F_of(shape(mid)) >= F:
            hi = mid
        else:
            lo = mid
    return hi


pitch, printed = 46.0, 75.8
N_fit = int(printed * 1000 / pitch)          # rows that fit if the rows had the whole printed length
res = dict(f_max=f_max, f_star_check=f_at(w_in), pitch_um=pitch, printed_mm=printed, rows_in_printed_length=N_fit, readings=[])
for Flab, F in (('outlet identity (the presets\' F)', 0.5592819159228324), ('assumed constant-shear taper reaching 52 um', 0.5105),
                ('measured (Dinh Fig. 4, mean)', 0.505), ('measured minus one SD', 0.472)):
    n_bound = rows_needed(F, lambda N: [w_min] * N)
    n_lin = rows_needed(F, linear_all)
    # the largest q (0..1, step 0.001) with linear_then_flat fitting in N_fit rows
    q_fit = None
    for k in range(1000, -1, -1):
        q = k / 1000.0
        if F_of(linear_then_flat(N_fit, q)) >= F:
            q_fit = q
            break
    mean_needed = (1.0 - (1.0 - F) ** (1.0 / N_fit)) / 2.0
    lo_w, hi_w = 20.0, w_in                  # the constant centre width that extracts mean_needed per row
    for _ in range(80):
        mid = 0.5 * (lo_w + hi_w)
        lo_w, hi_w = (mid, hi_w) if f_at(mid) > mean_needed else (lo_w, mid)
    row = dict(F_reading=Flab, F=round(F, 4), rows_every_row_at_52um=n_bound, length_mm_every_row_at_52um=round(n_bound * pitch / 1000, 2),
               rows_linear_75_to_52=n_lin, length_mm_linear_75_to_52=round(n_lin * pitch / 1000, 2),
               constant_f_needed_in_printed_length_over_f_max=round(mean_needed / f_max, 4),
               the_same_as_a_constant_centre_width_um=round(lo_w, 2),
               largest_share_of_rows_tapering_if_rows_had_the_whole_length=q_fit)
    res['readings'].append(row)
    print(json.dumps(row), flush=True)


def q_fit(F, n_rows, credit=0.0):
    for k in range(1000, -1, -1):
        q = k / 1000.0
        if F_of(linear_then_flat(n_rows, q), credit) >= F:
            return q
    return None


credits = {}
for k in ('depicted', 'specified'):
    jp = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dinh_length_bound_%s.json' % k)
    if os.path.exists(jp):
        credits[k] = float(json.load(open(jp))['rows_credited_to_this_head'])
res['sensitivity'] = []
for F in (0.5105, 0.505, 0.472):
    for pt in (45.0, 46.0, 47.0):
        n_fit = int(printed * 1000 / pt)
        row = dict(F=F, pitch_um=pt, rows_in_printed_length=n_fit, largest_share_tapering=q_fit(F, n_fit),
                   every_row_bound_mm=round(rows_needed(F, lambda N: [w_min] * N) * pt / 1000, 2))
        for k, c in credits.items():
            row['largest_share_tapering_with_%s_head_credit' % k] = q_fit(F, n_fit, c)
        res['sensitivity'].append(row)
        print('sensitivity', json.dumps(row), flush=True)
res['head_credit_rows'] = credits

# The cutoff law below its kappa plateau (after the second and third independent checks).  Everything above holds the
# lamina constant at every gap, with the plateau prefactor kappa_g: App. B's own reading of Dinh's rule ("Under the
# Xia-scaled reading, which is the one Dinh describes, the lamina is the same at every gap", App. B :2970-2972;
# glossary :1336-1342).  Under it kappa cancels from the row count; below the plateau (u = D / w_c < U_LO, i.e.
# w_c > 62.1 um at D 2.1742 um) only the cutoff's VALUE is uncertain: App. B's band for this reading is
# [1.80, 2.17] um ("the band collapses to", App. B :2972; its recorded P-C7.1 band, which holds f at its initial
# value instead, is [1.25, 2.17] um), the records' 1.805-2.174 um.  The readings that hold D
# itself constant and use the frozen Hill curve kappa(u) below the plateau are computed for comparison.  The Hill
# curve lies below kappa_g everywhere below the plateau (at most 4.418 at u = U_LO), so at a given D the lamina is
# thicker; but a D that returns f* at 75 um with it is smaller, and kappa(u) rises as w_c narrows, so the lamina thins
# downstream.  Those readings lengthen the device; the others do not return f* (after the fourth check):
#   * D at the lower end of that collapsed band, 1.805 um (the Hill kappa at the nominal u 0.029; it does not return
#     f* at 75 um): a longer every-row bound, a slightly shorter linear taper;
#   * D at the self-consistent Hill roots of D = kappa(D / 75) x(f*) 75, which do return f* at 75 um: longer;
#   * D 2.1742 um held with the Hill kappa, a reading a check proposed: it does not return f* at 75 um either; the
#     bound is unchanged and the linear taper shorter.
kmax, uh, nh = S.HILL


def kap(u):
    return S.KAPPA_G if u >= S.U_LO else kmax * u ** nh / (u ** nh + uh ** nh)


def g3(x):
    return 3 * x ** 2 - 2 * x ** 3


D_plateau = S.KAPPA_G * x_star * w_in
D_low = kap(D_plateau / w_in) * x_star * w_in               # App. B's collapsed band's lower end (not P-C7.1's 1.25)
hfun = lambda D: D - kap(D / w_in) * x_star * w_in
roots, prev = [], None
for i in range(50, 2200):
    Dv = i * 0.001
    v = hfun(Dv)
    if prev is not None and (v > 0) != (prev[1] > 0):
        lo_, hi_ = prev[0], Dv
        for _ in range(60):
            mid = 0.5 * (lo_ + hi_)
            if (hfun(mid) > 0) == (hfun(lo_) > 0):
                lo_ = mid
            else:
                hi_ = mid
        roots.append(0.5 * (lo_ + hi_))
    prev = (Dv, v)
readings = {'constant lamina, plateau prefactor at every gap (App. B\'s reading of Dinh\'s rule; the presets and everything above)': (D_plateau, lambda w: g3(D_plateau / (S.KAPPA_G * w)))}
readings['constant D at the lower end of App. B\'s band for a constant lamina (%.3f um; does not return the printed f* at 75 um), Hill kappa below the plateau' % D_low] = (D_low, lambda w: g3(D_low / (kap(D_low / w) * w)))
for r_ in roots:
    readings['constant D %.3f um, Hill kappa below the plateau (self-consistent: returns f* at 75 um)' % r_] = (r_, (lambda Dv: (lambda w: g3(Dv / (kap(Dv / w) * w))))(r_))
readings['constant D 2.1742 um held with the Hill kappa below the plateau (does not return the printed f* at 75 um)'] = (D_plateau, lambda w: g3(D_plateau / (kap(D_plateau / w) * w)))
res['cutoff_law_below_plateau'] = []
for lab, (Dv, ff) in readings.items():
    fm = ff(w_min)
    lnfm = math.log(1.0 - 2.0 * fm)

    def Fw(widths, ff=ff):
        return 1.0 - math.exp(sum(math.log(1.0 - 2.0 * ff(w)) for w in widths))

    def rows_n(F, shape, Fw=Fw):
        lo, hi = 2, 40000
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if Fw(shape(mid)) >= F:
                hi = mid
            else:
                lo = mid
        return hi

    def qf(F, n_rows, Fw=Fw):
        for k in range(1000, -1, -1):
            if Fw(linear_then_flat(n_rows, k / 1000.0)) >= F:
                return k / 1000.0
        return None
    row = dict(reading=lab, D_um=round(Dv, 4), f_at_75_um=ff(w_in), f_at_52_um=fm,
               every_row_bound_mm_at_F_0505=round(rows_n(0.505, lambda N: [w_min] * N) * pitch / 1000, 2),
               every_row_bound_mm_at_F_0472=round(rows_n(0.472, lambda N: [w_min] * N) * pitch / 1000, 2),
               linear_taper_mm_at_F_0505=round(rows_n(0.505, linear_all) * pitch / 1000, 2),
               largest_share_tapering_at_F_0505={'45_um': qf(0.505, int(printed * 1000 / 45.0)), '46_um': qf(0.505, N_fit),
                                                 '47_um': qf(0.505, int(printed * 1000 / 47.0))})
    res['cutoff_law_below_plateau'].append(row)
    print('below-plateau', json.dumps(row), flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
