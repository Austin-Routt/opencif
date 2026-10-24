"""The design code's reports and acceptance validators, for every opencif build.

The design code (reference/CIF_v2.py v2.1.0) attaches mandatory reports to every device its `design()` makes, and its
acceptance run (`run_acceptance`) records three operating validators. None of them is a pass condition there, and
none is one here. Each is computed on opencif's design (and, where it says so, on the drawn mask) by the design
code's own formula and is recorded in `validation.design_code_reports`. The reports and validators are as follows:

  B3  the pressure drop along the center channel, sum R_c(i) Q_c(i) at the design flow (design()'s reports.dP_tot),
      and the drawn device's pressure drop from gap 1 to the outlets, from the drawn network
  B6  the yield-corrected filtrate fraction F_real = eta F with eta = 0.903 (measured on Dinh's V1; design()'s F_real)
  B7  the transfer bias D_true/D_predicted = sqrt(mean c2(beta_cal) / c2(h/w_c(i))) over the rows
      (Pairing.transfer_bias, spec v5 S5)
  B8  the sensitivity map (spec v5 S16-e), which holds dln F/dln w_s(out), dln F/dln w_c(out), dln F/dln h,
      dln D/dln w_c, dln D/dln h, dln D/dln f, and dln N/dln kappa = 2, with the design code's own check M17
      (|dln F/dln h| < |dln F/dln w_s(out)|/3 when both outlet channels are narrower than h)
  B9  the wall porosity phi = G/Lambda (Process.phi), whereas the kernels assume a solid wall
  V_operating_regime  the gates Re <= 1, 0.05 Re D_h <= Lambda/2, the peak wall shear C(alpha) Q/(a b^2) <= 1500 1/s
      at the inlet, and the trunk's pressure drop sum r(w_c) Lambda P Q <= 200 kPa, at 5 uL/min and at the design flow
  V_steric  the limits D/G < 0.5, D/h < 0.3, and delta*/h < 0.2
  V_fab_tolerance  24 draws (seed 0) that add N(0, 0.25 um) to every center and side width, giving the share of draws
      whose F is within 3e-3 of the target (ok when >= 0.95) and the adherence (the share of the extracted flow at
      rows whose lamina is within 5% of delta*)

`tests/test_design_code.py` checks each formula against the design code's own function on the same inputs.
"""
import math

import numpy as np

from opencif.physics import kernels as K
from opencif.physics import network as NW
from opencif.physics import schedule as S

ETA = K.Pairing().eta                  # 0.903
Q_ACC_UL_MIN = 5.0                     # V_operating_regime's flow
RE_ACC, GDOT_ACC, DP_ACC_KPA = 1.0, 1500.0, 200.0
STERIC = dict(D_over_G=0.5, D_over_h=0.3, dstar_over_h=0.2)
FAB_N, FAB_SIGMA_UM, FAB_SEED, FAB_F_TOL, FAB_YIELD = 24, 0.25, 0, 3e-3, 0.95


def _f(x, n=4):
    """Round x to n significant figures for the record, or return None when x is None or not finite."""
    return None if x is None or not np.isfinite(x) else float('%.*g' % (n, x))


def sensitivity_map(wc_out_um, ws_out_um, wc_in_um, f_in, h_um, kappa, mu=S.MU):
    """Return the sensitivity map of spec v5 S16-e (the design code's sensitivity_map) as log-derivatives at the
    design point, taken by central differences with d = 1e-4.

    F follows the PE outlet identity F = 2/(r(w_s,out)/r(w_c,out) + 2), and D = kappa w_c x(f). The cubic's x(f)
    does not depend on h, so dln D/dln h = 0."""
    UM = S.UM
    h, wc_o, ws_o, wc_i = h_um * UM, wc_out_um * UM, ws_out_um * UM, wc_in_um * UM
    d = 1e-4

    def dln(fn, x):
        return (math.log(fn(x * (1 + d))) - math.log(fn(x * (1 - d)))) / (2 * d)

    F_ws = lambda ws: 2.0 / (K.r(ws, h, mu) / K.r(wc_o, h, mu) + 2.0)
    F_wc = lambda wc: 2.0 / (K.r(ws_o, h, mu) / K.r(wc, h, mu) + 2.0)
    F_h = lambda hh: 2.0 / (K.r(ws_o, hh, mu) / K.r(wc_o, hh, mu) + 2.0)
    D_wc = lambda wc: kappa * wc * K.x_of_f(f_in, h / wc)
    D_h_ = lambda hh: kappa * wc_i * K.x_of_f(f_in, hh / wc_i)
    D_f = lambda ff: kappa * wc_i * K.x_of_f(ff, h / wc_i)
    return {"dlnF/dln w_s(out)": dln(F_ws, ws_o), "dlnF/dln w_c(out)": dln(F_wc, wc_o), "dlnF/dln h": dln(F_h, h),
            "dlnD/dln w_c": dln(D_wc, wc_i), "dlnD/dln h": dln(D_h_, h), "dlnD/dln f": dln(D_f, f_in),
            "dlnN/dln kappa": 2.0}


def operating_regime(wc_rows_um, f_rows, h_um, pitch_um, q_ul_min, mu=1.0e-3, rho=1.0e3):
    """Evaluate the design code's V_operating_regime on a row schedule, with widths in micrometers. The validator
    checks the inlet's Re and naive entrance length, the peak wall shear at the inlet, and the trunk-only pressure
    drop at the flow."""
    UM = S.UM
    h, Lam, wc0 = h_um * UM, pitch_um * UM, wc_rows_um[0] * UM
    Q = q_ul_min * 1e-9 / 60.0
    Dh = 2.0 * wc0 * h / (wc0 + h)
    ub = Q / (wc0 * h)
    Re = rho * ub * Dh / mu
    Le = 0.05 * Re * Dh
    a0 = max(wc0, h) / min(wc0, h)
    gd = K.C_shear(a0) * Q / (max(wc0, h) * min(wc0, h) ** 2)
    P_after, dP = 1.0, 0.0
    for wc, f in zip(wc_rows_um, f_rows):
        dP += K.r(wc * UM, h, mu) * Lam * (P_after * Q)
        P_after *= (1.0 - 2.0 * f)
    g = dict(Re=bool(Re <= RE_ACC), entrance=bool(Le <= 0.5 * Lam), shear=bool(gd <= GDOT_ACC),
             pressure=bool(dP <= DP_ACC_KPA * 1e3))
    return dict(ok=all(g.values()), flow_uL_min=float(q_ul_min), Re=_f(Re), Le_um=_f(Le / UM), gdot_per_s=_f(gd),
                dP_kPa=_f(dP / 1e3), gates=g)


def steric(D_um, G_um, h_um, dstar_um):
    """Evaluate the design code's V_steric, the check of clogging and steric passage (current practice)."""
    r1, r2, r3 = D_um / G_um, D_um / h_um, dstar_um / h_um
    return dict(ok=bool(r1 < STERIC['D_over_G'] and r2 < STERIC['D_over_h'] and r3 < STERIC['dstar_over_h']),
                D_over_G=_f(r1), D_over_h=_f(r2), dstar_over_h=_f(r3), limits=STERIC)


def fab_tolerance(sc, info, n=FAB_N, sigma_um=FAB_SIGMA_UM, seed=FAB_SEED, F_tol=FAB_F_TOL):
    """Evaluate the design code's V_fab_tolerance on opencif's design network.

    Each of n draws moves every row's center width and side width by N(0, sigma), drawn in the design code's order
    (for each row, the center and then the side). It rescales each element's resistance to its new width (a head
    unit at its channel width, a U-turn at the mean of its two sides, and the outlet stretch at the last row's
    widths), leaves the gaps unchanged, and solves the network again. The yield is the share of draws with
    |F - F_target| <= F_tol, and the adherence of a draw is the share of its extracted flow at rows whose lamina
    x(f) w_c is within 5% of delta*."""
    spec, h = sc.spec, sc.spec.h
    N = sc.N
    rr = S.rr
    units = {u['k']: float(u['w']) for u in info['units']}
    loops = sc.info.get('loop_w', {}) or {}
    turns = set(sc.info.get('turns', []))
    dstar = float(spec.delta_star)
    rng = np.random.default_rng(seed)
    Fs, As = [], []
    for _ in range(n):
        dc, ds = np.zeros(N), np.zeros(N)
        for k in range(N):
            dc[k] = float(rng.normal(0.0, sigma_um))
            ds[k] = float(rng.normal(0.0, sigma_um))
        wc_new = np.maximum(np.asarray(sc.wc, float) + dc, 2.0)
        Rc = np.array([sc.Rc[k] * rr(0.5 * (wc_new[k] + wc_new[min(k + 1, N - 1)]), h) /
                       rr(0.5 * (sc.wc[k] + sc.wc[min(k + 1, N - 1)]), h) for k in range(N - 1)])
        Rs = np.array(sc.Rs, float).copy()
        for k in range(N - 1):
            if k in units:
                w0 = units[k]
            elif k in turns:
                w0 = float(loops.get(k, sc.ws[k - 1]))
            elif np.isfinite(sc.ws[k]):
                w0 = float(sc.ws[k])
            else:
                continue
            Rs[k] = sc.Rs[k] * rr(max(w0 + ds[k], 2.0), h) / rr(w0, h)
        wce, wse = float(sc.wc[-1]), float(sc.ws[N - 2])
        Rc_out = sc.Rc_out * rr(max(wce + dc[-1], 2.0), h) / rr(wce, h)
        Rs_out = sc.Rs_out * rr(max(wse + ds[-1], 2.0), h) / rr(wse, h)
        r = NW.solve(Rc, Rs, sc.Rg, Rc_out, Rs_out)
        Fs.append(float(r['F']))
        P, w_in, w_tot = 1.0, 0.0, 0.0
        for fe, wc in zip(np.asarray(r['f'], float), wc_new):
            fc = min(max(float(fe), 1e-9), 0.4999)
            dr = K.x_of_f(fc) * wc / dstar
            wgt = 2.0 * float(fe) * P
            w_tot += wgt
            if abs(dr - 1.0) <= 0.05:
                w_in += wgt
            P *= (1.0 - 2.0 * float(fe))
        As.append(w_in / w_tot if w_tot > 0 else 0.0)
    Fs, As = np.asarray(Fs), np.asarray(As)
    Fd = float(spec.F_target)
    yld = float((np.abs(Fs - Fd) <= F_tol).mean())
    return dict(ok=bool(yld >= FAB_YIELD), n=n, sigma_um=sigma_um, seed=seed, F_target=Fd, F_tol=F_tol,
                F_mean=_f(float(Fs.mean()), 6), F_std=_f(float(Fs.std()), 3), F_min=_f(float(Fs.min()), 6),
                F_max=_f(float(Fs.max()), 6), yield_F=yld, adherence_mean=_f(float(As.mean())),
                adherence_min=_f(float(As.min())))


def reports(sc, info, m, spec, F_drawn):
    """Compute every report and acceptance validator above for one build. Here sc is the design as drawn, m is the
    measurement of the drawn polygon, which gives the drawn device's pressure drop, and F_drawn enters B6."""
    from opencif.checks import measure as M
    h, G, g_p = float(spec.h), float(spec.G), float(spec.g_p)
    pitch = G + g_p
    wc, f = np.asarray(sc.wc, float), np.asarray(sc.f, float)
    D = S.KAPPA_G * float(spec.delta_star)
    Q = float(spec.Q_in_uL_min) * 1e-9 / 60.0
    mu, rho = float(getattr(spec, 'mu_f', S.MU)), float(getattr(spec, 'rho_f', 1000.0))
    Pin = np.concatenate([[1.0], np.asarray(sc.P, float)[:-1]])
    out = {}
    # B3 the pressure drop, from resistances at water's viscosity rescaled to the spec's fluid
    dp_trunk = (float(np.sum(np.asarray(sc.Rc, float)[:-1] * Pin[:-1] * Q)) + float(sc.Rc_out) * float(sc.P[-1]) * Q) \
        * mu / S.MU
    o = m['outlet']
    r = M.solve3(m['Rc'][:-1], m['RsA'], m['RsB'], m['RgA'], m['RgB'], o['Rc_out'], o['RsA_out'], o['RsB_out'])
    dp_drawn = float(r['Pc'][0]) * Q * mu / S.MU
    out['B3_pressure_drop'] = dict(
        trunk_kPa=_f(dp_trunk / 1e3), drawn_device_kPa=_f(dp_drawn / 1e3), flow_uL_min=float(spec.Q_in_uL_min),
        note='trunk: sum R_c(i) Q_c(i) over the design rows and the outlet stretch (the design code\'s dP_tot); '
             'drawn device: the inlet pressure of the drawn network from gap 1 to the outlets (the inlet stem, its '
             'filter and the off-chip tubing not counted); both at the spec\'s fluid')
    # B6 the design code's yield correction, which is reported and never applied
    out['B6_F_real'] = dict(eta=ETA, F_drawn=_f(F_drawn, 6), eta_F=_f(ETA * F_drawn, 6),
                            note='F_real = eta F, eta measured on Dinh\'s V1 (0.505 / 0.5592): reported, never applied '
                                 '(docs/LIMITATIONS.md 5)')
    # B7 the transfer bias
    lo_, hi_ = K.Pairing().transfer_bias([h / w for w in wc])
    out['B7_transfer_bias'] = dict(
        D_true_over_D_predicted=[_f(lo_), _f(hi_)], beta_cal=list(K.Pairing().beta_cal),
        note='sqrt(mean c2(beta_cal) / c2(h/w_c)) over the rows: the design code\'s mandatory report, never applied; '
             'beta_cal is the design code\'s assumed calibration aspect ratios ("ASSUMED (h_cal = 80 um).  '
             'UNVERIFIABLE."); App. B records the transfer bias as an uncertainty band rather than applying it '
             '(:3985-3986)')
    # B8 the sensitivity map, at the last row's widths (the outlet) and the first row's
    ws_out = float(sc.ws[sc.N - 2])
    sm = sensitivity_map(float(wc[-1]), ws_out, float(wc[0]), float(f[0]), h, S.KAPPA_G)
    m17 = (abs(sm['dlnF/dln h']) < abs(sm['dlnF/dln w_s(out)']) / 3.0) if (wc[-1] < h and ws_out < h) else None
    out['B8_sensitivity_map'] = dict(
        {k: _f(v) for k, v in sm.items()}, M17_widths_carry_the_device=m17,
        note='spec v5 S16-e at the outlet (the last centre and side widths) and at gap 1; M17 (the design code\'s '
             'check): |dlnF/dln h| < |dlnF/dln w_s(out)|/3 when both outlet channels are narrower than h (None: not '
             'applicable)')
    # B9 the wall porosity
    out['B9_wall_porosity'] = dict(phi=_f(G / pitch), note='G/Lambda: the kernels assume a solid wall (the design '
                                                            'code\'s Process.phi, reported)')
    # V_operating_regime, at the design code's 5 uL/min and at the design flow
    out['V_operating_regime'] = dict(at_5_uL_min=operating_regime(wc, f, h, pitch, Q_ACC_UL_MIN, mu, rho),
                                     at_design_flow=operating_regime(wc, f, h, pitch, float(spec.Q_in_uL_min), mu, rho),
                                     limits=dict(Re=RE_ACC, entrance='0.05 Re D_h <= Lambda/2', gdot_per_s=GDOT_ACC,
                                                 dP_kPa=DP_ACC_KPA))
    out['V_steric'] = steric(D, G, h, float(spec.delta_star))
    out['V_fab_tolerance'] = fab_tolerance(sc, info)
    out['role'] = ('reported, as in the design code: design()\'s mandatory reports and run_acceptance\'s operating '
                   'validators are not pass conditions there (docs/CIF_V2_REQUIREMENTS.md)')
    return out
