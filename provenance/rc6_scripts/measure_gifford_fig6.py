"""rc6: Gifford 2014 Fig. 6, the micrograph (right panel): the post geometry at its printed 125 um centre width.

Fig. 6's caption prints a 250 um scale bar, but the bar runs off the raster's right edge, so it bounds the scale from
above only (its visible length, the left end counted by coverage); the left panel is a schematic (topology only).
The caption's centre width w_c = 125 um is the one printed length the micrograph shows whole, so it is the scale, and
the bar is the check.  rc5 used Fig. 4's l_g 56.8 for this device and an unscripted pitch reading of 52.7 um; this script
reads both on Fig. 6 itself.

Method (all lengths in raster pixels until the last step; the raster is pdfimages' gif7-000.png of the project PDF):
  1. the panel (x 199-827, y 59-540) is resampled (bicubic, 0.25 px) into a frame with u along the post rows and v
     across them; the rows' direction is iterated until the apex lines of the posts are level (a fit of the apex
     positions against u), starting from the long straight walls' direction;
  2. the post rows are the v-bands where the darkness is periodic along u (the pills, 10-13 px apart);
  3. along each row, the pills' two side outlines (dark lines) give each pill's centre and width a (outline centre to
     outline centre); the pitch is the slope of the centres against their index (least squares);
  4. across each row, each pill's two apexes (its inner and outer tips) are located on the mean profile of the
     columns within +-1 px of its centre: the strongest darkness peak on each side of the row's middle, fitted with a
     Gaussian on a linear baseline (the outermost peak is not used: on the channel side the cells in the centre
     channel make weaker peaks beyond the tip); the apexes of a row lie on two tip lines (least squares, refitted
     without residuals beyond 3 robust SDs); l (apex to apex) is the post's length across the flow; the centre
     channel is the distance between the inner tip lines of the leg's two rows;
  5. synthetic calibration: the same steps run on a rendered image of the leg (stadium outlines at the measured
     geometry, drawn as Gaussian lines whose sigma is matched to the micrograph's straight walls through the same
     resampling and fit, integrated over the raster's pixels, at the same angle) give the estimators' biases, which
     are removed (the tips read ~0.09 px inside the outline's apex; pitch and width a are unbiased); the corrected
     geometry rendered and measured again returns the raw readings (closure), and the biases at 0.7x and 1.3x the
     matched sigma widen the ranges;
  6. scale = 125 um / the centre channel (each leg); pitch, l_g and the along-flow gap G_read = pitch - a in um.
The three legs with both post rows in the panel (the top leg's outer row is cut by the panel's edge) are measured
separately; the values are their mean and the ranges their extremes.  Checks reported, not used: the bar's bound;
the channel between the end blocks' straight edges (no tip estimator).

rc6 history: the first rc6 reading of this figure (2026-10-16, lost with its container before it was archived) gave
a pitch of 51.78 um and l_g 71.9 um.  Re-implementing it showed that taking the outermost peak as the tip, as it did,
picks a cell in the channel beside the lower rows' inner tips on part of the pills (their tip-line scatter 0.9 px
against 0.1-0.2 px on every other tip line), narrowing the channel by ~0.3 px and raising the scale by ~1 %; this
version takes the strongest peak.  g_p for the presets is the pitch minus G = 20.5 um (Fig. 4's printed gap; the micrograph's G_read is
reported, not used: relief contrast shifts outlines along the flow, as on Fig. 4's inset).
Usage: python3 -I measure_gifford_fig6.py GIF7_000_PNG OUT.json"""
import hashlib
import json
import sys

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy import optimize

RAS, OUT = sys.argv[1], sys.argv[2]
SHA = hashlib.sha256(open(RAS, 'rb').read()).hexdigest()
IM = np.asarray(Image.open(RAS).convert('L')).astype(float)
X0, X1, Y0, Y1 = 199, 827, 59, 540                     # the micrograph, one pixel inside its frame
BG = float(np.median(IM[Y0:Y1 + 1, X0:X1 + 1]))
W_C_UM = 125.0                                         # the printed centre width (Fig. 6 caption; p. 4502)
G_UM = 20.5                                            # the printed gap, Fig. 4 (the presets' G)
DU = 0.25


def frame(theta):
    th = np.radians(theta)
    eu = np.array([np.cos(th), -np.sin(th)])          # along the rows (up-right in the raster)
    ev = np.array([np.sin(th), np.cos(th)])           # across
    return eu, ev


def resample(img, theta, x0=X0, x1=X1, y0=Y0, y1=Y1, bg=None):
    eu, ev = frame(theta)
    C = np.array([[x0, y0], [x1, y0], [x0, y1], [x1, y1]], float)
    U, V = C @ eu, C @ ev
    us = np.arange(np.floor(U.min()), U.max(), DU)
    vs = np.arange(np.floor(V.min()), V.max(), DU)
    UU, VV = np.meshgrid(us, vs)
    X = UU * eu[0] + VV * ev[0]
    Y = UU * eu[1] + VV * ev[1]
    ok = (X >= x0 + 1) & (X <= x1 - 1) & (Y >= y0 + 1) & (Y <= y1 - 1)
    R = ndi.map_coordinates(img, [Y, X], order=3, mode='nearest')
    b = BG if bg is None else bg
    return np.where(ok, np.clip(b - R, 0.0, None), np.nan), us, vs


def periodic_power(sig, x, k_lo=1 / 13.0, k_hi=1 / 10.0, n=31):
    r = sig - sig.mean()
    ks = np.linspace(k_lo, k_hi, n)
    return max(abs(np.sum(r * np.exp(-2j * np.pi * k * x))) / len(r) for k in ks)


def rows_of(D, us, vs, min_u=150.0):
    """v-bands where the darkness is periodic along u: [(v_lo, v_hi)] of each post row"""
    P = np.zeros(len(vs))
    for i in range(len(vs)):
        ok = ~np.isnan(D[i])
        if ok.sum() * DU >= min_u:
            P[i] = periodic_power(D[i][ok], us[ok])
    on = P > 0.3 * P.max()
    bands, s = [], None
    for i, o in enumerate(np.append(on, False)):
        if o and s is None:
            s = i
        elif not o and s is not None:
            if (i - s) * DU >= 6.0:                       # a pill row is ~16 px across; edge slivers are dropped
                bands.append((vs[s], vs[i - 1]))
            s = None
    return bands


def gauss_peak(x, y, x0, half=1.25):
    m = np.abs(x - x0) <= half
    xx, yy = x[m], y[m]
    if len(xx) < 6 or np.any(np.isnan(yy)):
        return np.nan
    try:
        p, _ = optimize.curve_fit(lambda t, A, c, s, b0, b1: A * np.exp(-(t - c) ** 2 / (2 * s * s)) + b0 + b1 * (t - x0),
                                  xx, yy, p0=[yy.max() - yy.min(), xx[np.argmax(yy)], 0.8, yy.min(), 0.0], maxfev=4000)
    except RuntimeError:
        return np.nan
    return p[1] if abs(p[1] - x0) <= half and 0.2 < abs(p[2]) < 3.0 else np.nan


def peaks(x, y, min_sep=1.5):
    out = []
    for i in range(1, len(y) - 1):
        if y[i] >= y[i - 1] and y[i] > y[i + 1] and y[i] > 0.35 * np.nanmax(y):
            if not out or x[i] - out[-1] >= min_sep:
                out.append(x[i])
            elif y[i] > y[np.searchsorted(x, out[-1])]:
                out[-1] = x[i]
    return out


def robust_line(x, y, k=3.0):
    """least-squares line, refitted without points beyond k robust SDs (1.4826 MAD) of the residuals"""
    keep = np.ones(len(x), bool)
    for _ in range(5):
        f = np.polyfit(x[keep], y[keep], 1)
        r = y - np.polyval(f, x)
        sd = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        new = np.abs(r) <= max(k * sd, 0.05)
        if (new == keep).all():
            break
        keep = new
    return f, keep


def analyse_row(D, us, vs, band):
    """pills of one row: centres, widths a, pitch, apex lines (v_in side is decided by the caller)"""
    v_lo, v_hi = band
    iv = (vs >= v_lo - 6) & (vs <= v_hi + 6)               # the band is the periodic core; the tips lie ~2-3 px beyond
    vmid0, vmid1 = v_lo + 0.3 * (v_hi - v_lo), v_hi - 0.3 * (v_hi - v_lo)
    im_ = (vs >= vmid0) & (vs <= vmid1)
    s = np.nanmean(D[im_], axis=0)
    ok = ~np.isnan(s)
    # the posts' u-extent: local periodic power over 3 pitches
    w = int(round(36 / DU))
    jj = np.arange(w, len(us) - w, 4)
    pj = np.array([0.0 if np.any(np.isnan(s[j - w:j + w])) else periodic_power(s[j - w:j + w], us[j - w:j + w]) for j in jj])
    pw = np.interp(np.arange(len(us)), jj, pj, left=0.0, right=0.0)
    on = pw > 0.4 * pw.max()
    runs, st = [], None
    for j, o in enumerate(np.append(on, False)):
        if o and st is None:
            st = j
        elif not o and st is not None:
            runs.append((st, j - 1))
            st = None
    j0, j1 = max(runs, key=lambda r: r[1] - r[0])
    u_lo, u_hi = us[j0] - w * DU + 6.0, us[j1] + w * DU - 6.0     # the run's ends, less half a pitch
    u_lo, u_hi = max(u_lo, us[ok].min() + 2.0), min(u_hi, us[ok].max() - 2.0)
    m = (us >= u_lo) & (us <= u_hi) & ok
    pk = peaks(us[m], s[m])
    pk = [p for p in (gauss_peak(us[m], s[m], p0) for p0 in pk) if not np.isnan(p)]
    # pair the side outlines: within a pill the separation is the larger one (a > pitch - a here)
    d = np.diff(pk)
    big = d > np.median(d)
    cent, wid = [], []
    for i in range(len(d)):
        if big[i]:
            cent.append(0.5 * (pk[i] + pk[i + 1]))
            wid.append(d[i])
    cent, wid = np.array(cent), np.array(wid)
    # index the pills by their pitch (a missing pill keeps its index)
    p0 = np.median(np.diff(cent))
    idx = np.round((cent - cent[0]) / p0).astype(int)
    fit = np.polyfit(idx, cent, 1)
    pitch, res = fit[0], cent - np.polyval(fit, idx)
    # apexes: each pill's +-1 px columns, the profile across v; its two extreme peaks
    apex = []
    for uc in cent:
        cm = (us >= uc - 1.0) & (us <= uc + 1.0)
        prof = np.nanmean(D[np.ix_(iv, cm)], axis=1)
        vv = vs[iv]
        # the tips: the strongest darkness peak in each half of the window (the centre channel's cells make weaker
        # peaks beside the inner tip, so the outermost peak is not used)
        vc_ = 0.5 * (v_lo + v_hi)
        lo_h, hi_h = vv < vc_, vv > vc_
        if not lo_h.any() or not hi_h.any() or np.all(np.isnan(prof)):
            apex.append((uc, np.nan, np.nan))
            continue
        lo_, hi_ = (gauss_peak(vv, prof, vv[h][int(np.nanargmax(prof[h]))]) for h in (lo_h, hi_h))
        apex.append((uc, lo_, hi_))
    apex = np.array(apex)
    good = ~np.isnan(apex[:, 1]) & ~np.isnan(apex[:, 2])
    if good.sum() < 6:
        return None
    fl, gl = robust_line(apex[good, 0], apex[good, 1])
    fh, gh = robust_line(apex[good, 0], apex[good, 2])
    return dict(band=band, u_range=(float(u_lo), float(u_hi)), n_pills=int(len(cent)), pitch_px=float(pitch),
                pitch_resid_rms_px=float(res.std()), a_px=float(np.mean(wid)), a_sd_px=float(np.std(wid)),
                apex_lo=fl.tolist(), apex_hi=fh.tolist(), n_apex=int(good.sum()), n_apex_kept=[int(gl.sum()), int(gh.sum())],
                l_px=float(np.mean((apex[good, 2] - apex[good, 1])[gl & gh])),
                apex_lo_rms=float(np.std((apex[good, 1] - np.polyval(fl, apex[good, 0]))[gl])),
                apex_hi_rms=float(np.std((apex[good, 2] - np.polyval(fh, apex[good, 0]))[gh])))


def legs_of(rows):
    """pair consecutive rows whose facing tip lines are a centre channel apart (the narrow gap, ~27 px)"""
    legs, i = [], 0
    while i < len(rows) - 1:
        a, b = rows[i], rows[i + 1]
        gap = b['band'][0] - a['band'][1]
        if 15.0 < gap < 40.0:
            legs.append((a, b))
            i += 2
        else:
            i += 1
    return legs


def leg_values(a, b):
    """the centre channel between row a's upper (hi) tip line and row b's lower (lo) tip line, at the rows' common u"""
    u0, u1 = max(a['u_range'][0], b['u_range'][0]), min(a['u_range'][1], b['u_range'][1])
    um = 0.5 * (u0 + u1)
    wc = np.polyval(b['apex_lo'], um) - np.polyval(a['apex_hi'], um)
    return dict(u_mid=float(um), wc_px=float(wc), l_px=[a['l_px'], b['l_px']], a_px=[a['a_px'], b['a_px']],
                pitch_px=[a['pitch_px'], b['pitch_px']], slopes=[a['apex_lo'][0], a['apex_hi'][0], b['apex_lo'][0], b['apex_hi'][0]])


def measure(img, theta, bg=None):
    D, us, vs = resample(img, theta, bg=bg)
    rows = [analyse_row(D, us, vs, bd) for bd in rows_of(D, us, vs)]
    rows = [r for r in rows if r is not None and r['n_pills'] >= 8 and r['n_apex'] >= 6]
    return rows, legs_of(rows)


def wall_sigma(D, us, vs, v_guess):
    """the median Gaussian sigma (px, across v) of a straight wall's profile over 20 px chunks along u"""
    sig = []
    for u0 in np.arange(us.min(), us.max(), 20.0):
        cm = (us >= u0) & (us < u0 + 20.0)
        vm = np.abs(vs - v_guess) <= 4.0
        blk = D[np.ix_(vm, cm)]
        if np.isnan(blk).any() or blk.size == 0:
            continue
        prof = blk.mean(axis=1)
        vv = vs[vm]
        try:
            p, _ = optimize.curve_fit(lambda t, A, c, s_, b: A * np.exp(-(t - c) ** 2 / (2 * s_ * s_)) + b, vv, prof,
                                      p0=[prof.max(), vv[np.argmax(prof)], 0.8, 0.0], maxfev=4000)
            if abs(p[1] - v_guess) < 3.0 and 0.2 < abs(p[2]) < 3.0:
                sig.append(abs(p[2]))
        except RuntimeError:
            pass
    return float(np.median(sig)) if sig else np.nan, len(sig)


YY, XX = np.mgrid[0:IM.shape[0], 0:IM.shape[1]]
SS = 6                                                 # supersampling per raster pixel (box integration)
OFFS = (np.arange(SS) + 0.5) / SS - 0.5


def render(dist_fn, A, sig_line):
    """darkness A exp(-d^2 / 2 sigma^2) of the distance d to the outline, integrated over each raster pixel"""
    acc = np.zeros(IM.shape)
    for oy in OFFS:
        for ox in OFFS:
            d = dist_fn(XX + ox, YY + oy)
            acc += A * np.exp(-d * d / (2.0 * sig_line * sig_line))
    return BG - acc / (SS * SS)


def stadium_dist(u_rel, v_rel, a, l):
    """distance to the outline of a stadium (width a along u, apex-to-apex l along v), centred at the origin"""
    r = a / 2.0
    h = max(l / 2.0 - r, 0.0)
    vc = np.clip(v_rel, -h, h)
    return np.abs(np.hypot(u_rel, v_rel - vc) - r)


def leg_image(theta, a_px, l_px, pitch_px, v_c_rows, u_phase, sig_line, A):
    eu, ev = frame(theta)

    def dist(x, y):
        u = x * eu[0] + y * eu[1]
        v = x * ev[0] + y * ev[1]
        k = np.round((u - u_phase) / pitch_px)
        ur = u - (u_phase + k * pitch_px)
        return np.minimum.reduce([stadium_dist(ur, v - vc, a_px, l_px) for vc in v_c_rows])
    return render(dist, A, sig_line)


def synthetic_estimate(theta, geom, place, sig_line, A):
    """the estimators on a rendered leg of geometry geom = (a, l, pitch, centre) in px; place = (v of the upper row's
    inner tip line, band of the upper row) of the real leg"""
    a_px, l_px, p_px, wc = geom
    vin_a = place[0]
    img = leg_image(theta, a_px, l_px, p_px, (vin_a - l_px / 2.0, vin_a + wc + l_px / 2.0), 0.31, sig_line, A)
    rows_s, legs_s = measure(img, theta, bg=BG)
    best = min(legs_s, key=lambda q: abs(q[0]['band'][1] - place[1]))
    Ls = leg_values(*best)
    return dict(wc=Ls['wc_px'], l=float(np.mean(Ls['l_px'])), a=float(np.mean(Ls['a_px'])), pitch=float(np.mean(Ls['pitch_px'])))


def bar_lower_bound():
    """Fig. 6's 250 um bar (caption) runs off the raster's right edge: its visible length, the left end counted by
    coverage, is a lower bound on its length in px, so 250 um over it bounds the scale from above"""
    rows_ = [y for y in range(541, IM.shape[0]) if (IM[y, 780:IM.shape[1]] < 30).all()]
    blk = IM[min(rows_):max(rows_) + 1, :]
    prof = blk.mean(axis=0)
    dark = float(np.median(prof[-20:]))
    white = 255.0
    cov = np.clip((white - prof) / (white - dark), 0.0, 1.0)
    full = np.nonzero(cov > 0.98)[0]
    x_first = int(full.min())
    right_edge_dark = bool(cov[-1] > 0.98)
    length_min = float(cov[x_first - 1] + (IM.shape[1] - x_first))
    return dict(rows=[int(min(rows_)), int(max(rows_))], left_end_px=x_first - 1, left_end_coverage=round(float(cov[x_first - 1]), 3),
                runs_off_the_right_edge=right_edge_dark, visible_length_px=round(length_min, 3),
                um_per_px_upper_bound=round(250.0 / length_min, 4))


def main():
    # 1. the rows' direction: from the straight walls' 49.8 deg (a line trace), iterated until the tip lines are level
    theta = 49.8
    for _ in range(4):
        rows, legs = measure(IM, theta)
        sl = np.median([q for a_, b_ in legs for q in leg_values(a_, b_)['slopes']])
        theta_new = theta - np.degrees(np.arctan(sl))
        done = abs(theta_new - theta) < 0.002
        theta = theta_new
        if done:
            break
    rows, legs = measure(IM, theta)
    LEGS = [leg_values(a_, b_) for a_, b_ in legs]
    # 2. the line sigma, matched on the straight walls through the same resampling and fit
    D0, us0, vs0 = resample(IM, theta)
    mrow = np.nanmean(D0, axis=1)
    cnt = np.sum(~np.isnan(D0), axis=1) * DU
    cand = [i for i in range(3, len(vs0) - 3) if cnt[i] > 300 and mrow[i] == np.nanmax(mrow[i - 3:i + 4])
            and mrow[i] > 40.0 and periodic_power(D0[i][~np.isnan(D0[i])], us0[~np.isnan(D0[i])]) < 2.0]
    walls = [(float(vs0[i]),) + wall_sigma(D0, us0, vs0, vs0[i]) for i in cand]
    walls = [w for w in walls if w[2] >= 8]
    sigma_obs = float(np.median([w[1] for w in walls]))
    A = 1.6 * float(np.nanmax(mrow))
    ev = frame(theta)[1]
    v0 = float(np.mean([w[0] for w in walls])) + 0.37                   # an off-grid position

    def sigma_syn(sig_line):
        img = render(lambda x, y: np.abs(x * ev[0] + y * ev[1] - v0), A, sig_line)
        D_, u_, v_ = resample(img, theta, bg=BG)
        return wall_sigma(D_, u_, v_, v0)[0]
    sig_line = optimize.brentq(lambda q: sigma_syn(q) - sigma_obs, 0.2, 1.6, xtol=1e-3)
    # 3. synthetic calibration per leg (bias = estimate - truth at the measured geometry), its sensitivity to the line
    #    sigma (x0.7, x1.3), and a closure check: the corrected geometry rendered and measured returns the raw values
    per_leg = []
    for (ra, rb), L in zip(legs, LEGS):
        raw = dict(wc=L['wc_px'], l=float(np.mean(L['l_px'])), a=float(np.mean(L['a_px'])), pitch=float(np.mean(L['pitch_px'])))
        place = (float(np.polyval(ra['apex_hi'], L['u_mid'])), ra['band'][1])
        geom0 = (raw['a'], raw['l'], raw['pitch'], raw['wc'])
        bias = {}
        for f_ in (1.0, 0.7, 1.3):
            est = synthetic_estimate(theta, geom0, place, sig_line * f_, A)
            bias[f_] = {k: est[k] - v for k, v in zip(('a', 'l', 'pitch', 'wc'), geom0)}
        b = bias[1.0]
        corr = {k: raw[k] - b[k] for k in raw}
        clo = synthetic_estimate(theta, (corr['a'], corr['l'], corr['pitch'], corr['wc']), place, sig_line, A)
        closure = {k: clo[k] - raw[k] for k in raw}

        def um(c):
            s_ = W_C_UM / c['wc']
            return dict(um_per_px=s_, row_pitch_um=s_ * c['pitch'], l_g_um=s_ * c['l'], a_um=s_ * c['a'],
                        G_read_um=s_ * (c['pitch'] - c['a']))
        v1 = um(corr)
        alt = [um({k: raw[k] - bias[f_][k] for k in raw}) for f_ in (0.7, 1.3)]
        per_leg.append(dict(raw_px=raw, bias_px=b, corrected_px=corr, closure_px=closure, **v1,
                            sigma_sensitivity={k: [min(q[k] for q in alt), max(q[k] for q in alt)] for k in v1},
                            rows=[ra['band'], rb['band']], n_pills=[ra['n_pills'], rb['n_pills']]))

    # a cross-check that uses no tip estimator: each leg's centre channel continues between the two rows' end blocks,
    # whose facing edges are straight lines (a symmetric line profile, so a Gaussian fit is unbiased); the blocks need
    # not be flush with the post tips, so this is reported, not used
    D1, us1, vs1 = resample(IM, theta)
    blocks = []
    for (ra, rb), L in zip(legs, LEGS):
        u0 = max(ra['u_range'][0], rb['u_range'][0]) - 55.0
        u1 = min(ra['u_range'][0], rb['u_range'][0]) - 10.0
        cm = (us1 >= u0) & (us1 <= u1)
        e = []
        for vg in (np.polyval(ra['apex_hi'], L['u_mid']), np.polyval(rb['apex_lo'], L['u_mid'])):
            vm = np.abs(vs1 - vg) <= 4.0
            prof = np.nanmean(D1[np.ix_(vm, cm)], axis=1)
            vv = vs1[vm]
            e.append(gauss_peak(vv, prof, vv[int(np.nanargmax(prof))], half=1.5))
        blocks.append(dict(u_range=[float(u0), float(u1)], edges_v=[float(x) for x in e], separation_px=float(e[1] - e[0]),
                           um_per_px_if_125=float(W_C_UM / (e[1] - e[0]))))

    def mean_rng(k):
        v = [q[k] for q in per_leg]
        lo = min(min(v), min(q['sigma_sensitivity'][k][0] for q in per_leg))
        hi = max(max(v), max(q['sigma_sensitivity'][k][1] for q in per_leg))
        return float(np.mean(v)), [float(lo), float(hi)]
    s_m, s_r = mean_rng('um_per_px')
    p_m, p_r = mean_rng('row_pitch_um')
    l_m, l_r = mean_rng('l_g_um')
    g_m, g_r = mean_rng('G_read_um')
    a_m, a_r = mean_rng('a_um')
    pitch_px = float(np.mean([q['corrected_px']['pitch'] for q in per_leg]))
    bar = bar_lower_bound()
    out = dict(figure='Gifford 2014 Fig. 6, the micrograph (right panel)', raster='gif7-000.png', raster_sha256=SHA,
               panel_px=[X0, Y0, X1, Y1],
               scale_from='the printed centre width w_c = 125 um (Fig. 6 caption; p. 4502), between the inner tip lines '
                          'of each leg\'s two post rows',
               scale_bar=dict(bar, printed_um=250.0, note='the caption\'s 250 um bar runs off the raster\'s right edge; '
                              'its visible length bounds the scale from above only'),
               rows_direction_deg=round(theta, 3),
               line_sigma_px=dict(observed_on_walls=round(sigma_obs, 3), matched_line_sigma=round(sig_line, 3),
                                  walls=[[round(w[0], 2), round(w[1], 3), w[2]] for w in walls]),
               um_per_px=round(s_m, 4), um_per_px_range=[round(v, 4) for v in s_r],
               pitch_px=round(pitch_px, 4), row_pitch_um=round(p_m, 3), row_pitch_range_um=[round(v, 2) for v in p_r],
               row_pitch_upper_bound_um=round(pitch_px * bar['um_per_px_upper_bound'], 2),
               l_g_um=round(l_m, 2), l_g_range_um=[round(v, 1) for v in l_r],
               post_along_flow_um=round(a_m, 2), post_along_flow_range_um=[round(v, 1) for v in a_r],
               G_read_um=round(g_m, 2), G_read_range_um=[round(v, 1) for v in g_r],
               G_assumed_um=G_UM, g_p_um=round(p_m - G_UM, 2),
               preset_values=dict(post_length_um=round(p_m - G_UM, 1), post_width_um=round(l_m, 1)),
               ranges='the legs\' extremes, widened by the calibration at 0.7x and 1.3x the matched line sigma',
               end_blocks_check=dict(per_leg=blocks, note='the channel between each leg\'s two end blocks (straight '
                                     'edges, no tip estimator): %.2f-%.2f px, against %.2f-%.2f px between the corrected '
                                     'tip lines; reported, not used (the blocks need not be flush with the tips)'
                                     % (min(b_['separation_px'] for b_ in blocks), max(b_['separation_px'] for b_ in blocks),
                                        min(q['corrected_px']['wc'] for q in per_leg), max(q['corrected_px']['wc'] for q in per_leg))),
               per_leg=per_leg, rows=rows,
               note='g_p = pitch - G with G 20.5 um as Fig. 4 (printed); the micrograph\'s own G_read is reported, not '
                    'used: relief contrast shifts the outlines along the flow (as on Fig. 4\'s inset).  The left panel '
                    'is a schematic and is not used.')
    json.dump(out, open(OUT, 'w'), indent=1, default=float)
    print('rows direction %.3f deg; line sigma observed %.3f px, matched %.3f px; %d legs; bar >= %.3f px -> <= %.4f um/px'
          % (theta, sigma_obs, sig_line, len(per_leg), bar['visible_length_px'], bar['um_per_px_upper_bound']))
    for i, q in enumerate(per_leg):
        print('leg %d: centre %.3f px (raw %.3f, bias %+.3f, closure %+.4f) -> %.4f um/px; pitch %.4f px = %.2f um; l_g %.2f '
              'um; a %.2f um; G_read %.2f um' % (i + 1, q['corrected_px']['wc'], q['raw_px']['wc'], q['bias_px']['wc'],
                                                 q['closure_px']['wc'], q['um_per_px'], q['corrected_px']['pitch'],
                                                 q['row_pitch_um'], q['l_g_um'], q['a_um'], q['G_read_um']))
    print('end blocks: ' + out['end_blocks_check']['note'])
    print('scale %.4f um/px (%.4f-%.4f); pitch %.2f um (%.2f-%.2f; <= %.2f by the bar); l_g %.2f um (%.1f-%.1f); G_read %.2f; '
          'g_p = pitch - 20.5 = %.2f' % (s_m, s_r[0], s_r[1], p_m, p_r[0], p_r[1], out['row_pitch_upper_bound_um'], l_m,
                                         l_r[0], l_r[1], g_m, p_m - G_UM))


if __name__ == '__main__':
    main()
