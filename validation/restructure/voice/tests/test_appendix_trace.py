"""test_appendix_trace.py: each result of Appendices A and B, and each glossary value, that the generator relies on,
pinned to the code that implements it.

The suite is fast because it draws nothing. Run it with `python3 tests/test_appendix_trace.py`. Each check names the
result it pins by its App. A Result and equation number (as printed in APPENDIX_A_revised.pdf), by its App. B claim and
premise, or by its glossary macro. A failure means that the code no longer implements the frozen statement.
docs/ALIGNMENT.md maps every statement to its code, its check, and its status."""
import copy, math, os, re, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from opencif import spec as SP, build as BU, pipeline as PL   # noqa: F401
from opencif import pipeline as _PL_CODE                # noqa: E402  (every suite's log names the code it ran on)
print('opencif %s (code %s)' % (_PL_CODE.CODE['version'], _PL_CODE.CODE['sha256']))
from opencif.physics import schedule as S
from opencif.checks import measure as M2
from opencif.checks import design_code as DCR
sys.path.insert(0, os.path.join(ROOT, 'reference'))
import CIF_v2 as C                      # noqa: E402  the frozen design code

fails = []
def check(name, cond, info=''):
    print('%-78s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)

import json
GLOSSARY = json.load(open(os.path.join(ROOT, 'reference', 'glossary_values.json')))
PRE_PATH = os.path.join(ROOT, 'reference_glossary', 'PREAMBLE_FINAL.tex')
if os.path.exists(PRE_PATH):            # the dissertation's overlay, which holds the frozen preamble itself
    pre = open(PRE_PATH).read()
    mac = {m.group(1): m.group(2) for m in re.finditer(r'\\newcommand\{\\(\w+)\}\{([0-9.]+)\}', pre)}
else:                                   # the public repository, which holds the values that the checks below use
    pre = None
    mac = dict(GLOSSARY['macros'])

# ---- App. A Result 1 (eq. 10): the skim-flux cubic f(x) = 3x^2 - 2x^3 with x = delta/w_c -----------------------------
xs = np.linspace(0.0, 0.499, 51)
check('A R1 (eq. 10): f(x) = 3x^2 - 2x^3 (CIF_v2.g_cubic)', max(abs(C.g_cubic(x) - (3 * x * x - 2 * x ** 3)) for x in xs) < 1e-15)
check('A R1 (eq. 10): the design kernel is the cubic (schedule.f_of)',
      max(abs(S.f_of(x * 100.0, 100.0) - (3 * x * x - 2 * x ** 3)) for x in xs) < 1e-15)
# Step 3 of the section "For implementers" in App. A says to stop with an error if delta*/w_c >= 1/2.
try:
    S.f_of(50.0, 100.0); ok = False
except ValueError:
    ok = True
check('A implementers step 3: delta*/w_c >= 1/2 is an error (schedule.f_of)', ok)
# ---- App. A Result 2 (eqs. 18 and 20): the inversion x(f) ------------------------------------------------------------
fs = [1e-6, 1e-4, 3.34e-4, 1e-2, 0.1, 0.4]
check('A R2 (eq. 20): x(f) is the exact root (round trip)', max(abs(C.g_cubic(S.x_of(f)) - f) / f for f in fs) < 1e-12)
check('A R2 (eq. 18): leading order x = sqrt(f/3) for small f', abs(S.x_of(1e-6) / math.sqrt(1e-6 / 3) - 1) < 1e-3)
tv = ((1e-4, 0.005784668), (5.8e-4, 0.013969638), (0.1, 0.195800106))
dv_ = max(abs(S.x_of(f) - x) for f, x in tv)
check('A implementers test values: exact roots x(1e-4), x(5.8e-4), x(0.1)', dv_ < 6e-10, '%.1e' % dv_)
# The exact root is compared with the leading order over kappa_g's calibration range (the docs' 0.47% statement).
lo = max(abs(S.x_of(f) / math.sqrt(f / 3) - 1) for f in np.linspace(1e-5, 5.76e-4, 50))
check('A R2: exact root vs leading order <= 0.47 % up to the calibration cap 5.76e-4', lo < 0.0047, '%.3f %%' % (100 * lo))
# ---- App. A Results 5 and 6 (eq. 36) and the glossary pairing rule: D = kappa_g delta* with the cubic pair -----------
P = C.Pairing()
check('glossary: kappa_g = 4.9138, 95 %% CI [4.750, 5.078] (half-width %s %%)' % mac.get('kappCIhalfPct'),
      P.kappa == 4.9138 and tuple(P.kappa_ci) == (4.750, 5.078))
check('glossary pairing rule (A R9 corollary, eq. 59): kappa_g with the cubic, never with F_lam',
      P.kernel == 'cubic' and S.KAPPA_G == P.kappa)
kp = BU.KERNEL_PAIR
import inspect
src_b, src_p, src_a = inspect.getsource(BU.build), inspect.getsource(PL._write_refusal), open(os.path.join(ROOT, 'opencif', 'api.py')).read()
check('App. B C8 kernel declaration: the pair (build.KERNEL_PAIR) is the code\'s, and build, refusal and check records carry it',
      kp['kernel'] == 'cubic' and kp['kappa_g'] == S.KAPPA_G and kp['kappa_ci'] == list(S.KAPPA_CI)
      and 'kernel_pair=KERNEL_PAIR' in src_b and 'kernel_pair=BU.KERNEL_PAIR' in src_p and 'kernel_pair=BU.KERNEL_PAIR' in src_a)
# ---- The validity envelope -------------------------------------------------------------------------------------------
check('glossary envelope / App. B C8 CQ7: plateau u_lo 0.035, foot u_ff 0.012 (external inputs)',
      S.U_LO == 0.035 and S.U_FF == 0.012 and P.u_lo == 0.035 and P.u_ff == 0.012)
T = C.Tolerances()
check('App. B thesis scope / CIF_v2 M13, M14: gap dominance 0.05, Re_max 10, Re_stokes 1e-2 (build)',
      BU.M13_MAX == T.gap_dom == 0.05 and BU.RE_MAX == T.Re_max == 10.0 and BU.RE_STOKES == T.Re_stokes == 1e-2)
# ---- App. A Result 14 (eq. 88; Gifford 2014 eqn 5): the resistance of a rectangular duct -----------------------------
def r88(w, d, mu=S.MU, n=401):
    a, b = max(w, d), min(w, d)
    ssum = sum(math.tanh(k * math.pi * a / (2 * b)) / k ** 5 for k in range(1, n, 2))
    return 12 * mu / (a * b ** 3) / (1 - 192 * b / (math.pi ** 5 * a) * ssum)
devs = [abs(S.rr(w, h) / r88(w * S.UM, h * S.UM) - 1) for w, h in ((19, 100), (28, 140), (100, 80), (600, 140), (20, 20))]
check('A R14 (eq. 88): duct resistance = the Poiseuille series (5 aspect ratios)', max(devs) < 1e-6, '%.1e' % max(devs))
sq = S.rr(20.0, 20.0) * (20.0 * S.UM) ** 4 / S.MU
check('A implementers / glossary \\squareDuctR: square duct R h^4/(mu L) = 28.454154 (%.6f)' % sq, abs(sq - 28.454154) < 5e-7
      and abs(sq - float(mac['squareDuctR'])) < 5e-5)
# ---- App. A Result 10 (eqs. 60 to 67): the march and the implementers' march example ---------------------------------
# The example holds f = 1e-4 at every row of a constant-width channel and reads P, sigma, and F after 5,000 rows.
dm = SP.load(os.path.join(ROOT, 'examples', 'gifford_physics_linear.yaml'))
wc = dm['centre_channel']['width_in_um']
dm['targets']['cutoff_um'] = S.KAPPA_G * S.x_of(1e-4) * wc
dm['targets']['filtrate_fraction'] = 0.64
scm = S.plan(SP.to_engine(dm)[0])
P5, S5 = scm.P[4999], scm.S[4999]
sig = S5 / P5
check('A implementers march example: f 1e-4, 5000 rows: P 0.367843, sigma 0.859277, F 0.632157',
      abs(scm.f[0] - 1e-4) < 1e-15 and round(P5, 6) == 0.367843 and round(sig, 6) == 0.859277 and round(1 - P5, 6) == 0.632157
      and round(2 * sig / (1 + 2 * sig), 6) == 0.632157, 'P %.6f sigma %.6f' % (P5, sig))
# ---- The designs used below ------------------------------------------------------------------------------------------
base = SP.load(os.path.join(ROOT, 'examples', 'gifford_physics_linear.yaml'))
sc1 = S.plan(SP.to_engine(base)[0])
d2 = copy.deepcopy(base); d2['targets']['cutoff_um'] = 4.0
sc2 = S.plan(SP.to_engine(d2)[0])
for name_, sc_ in (('gifford_physics_linear', sc1), ('the same at D 4 um', sc2)):
    Pa, Sa = np.asarray(sc_.P), np.asarray(sc_.S)
    mi = np.max(np.abs(Pa + 2 * Sa - 1 - 2 * float(getattr(sc_.spec, 'S0', 0.0) or 0.0)))
    check('A R10 mass identity: |P + 2S - (1 + 2 sigma_0)| at every row (%s)' % name_, mi < 1e-12, '%.1e' % mi)
    phi = C.critical_phi(list(sc_.f[:sc_.N]))
    tl = abs(2 * phi - (1 - Pa[sc_.N - 1]))
    check('A R12 (eq. 80) / App. B C9 transport lemma: 2 phi_crit = F (%s)' % name_, tl < 1e-6, '%.1e' % tl)
    lam = np.asarray(sc_.delta[:sc_.N]) / sc_.spec.delta_star - 1
    check('A R13(a) (eq. 81) / App. B P-C8.1: the lamina at every row is delta* (%s)' % name_, np.max(np.abs(lam)) < 1e-9,
          '%.1e' % np.max(np.abs(lam)))
# ---- App. A Result 13(c) (eq. 86): D and F are independent inputs (forward substitution) -----------------------------
check('A R13(c) (eq. 86): changing D (hence f_gap) changes N, not F', abs(sc1.F - sc2.F) < 2 * max(np.max(sc1.f), np.max(sc2.f))
      and sc1.N != sc2.N, 'F %.5f / %.5f, N %d / %d' % (sc1.F, sc2.F, sc1.N, sc2.N))
d4 = copy.deepcopy(base); d4['targets']['filtrate_fraction'] = 0.5
sc4 = S.plan(SP.to_engine(d4)[0])
cut4 = S.KAPPA_G * np.array([S.x_of(f) for f in sc4.f[:sc4.N - 1]]) * np.asarray(sc4.wc[:sc4.N - 1])
check('A R13(c) (eq. 86): changing F changes N, not the cutoff', sc4.N < sc1.N and
      np.max(np.abs(cut4 / (S.KAPPA_G * sc4.spec.delta_star) - 1)) < 1e-9, 'N %d -> %d' % (sc1.N, sc4.N))
# ---- App. A Result 13(b) (eqs. 82 and 84) and App. B C8 P-C8.2: N D^2 is invariant at fixed F ------------------------
ND2 = []
for D in (3.0, 4.0, 5.0, 6.0, 7.0):
    dd = copy.deepcopy(base); dd['targets']['cutoff_um'] = D
    s_ = S.plan(SP.to_engine(dd)[0])
    ND2.append(s_.N * D * D)
F_ = base['targets']['filtrate_fraction']
closed = -S.KAPPA_G ** 2 * wc ** 2 * math.log(1 - F_) / 6.0          # the constant-w_c limit of eq. 82
spread = (max(ND2) - min(ND2)) / np.mean(ND2)
off = max(abs(v / closed - 1) for v in ND2)
check('A R13(b) / App. B P-C8.2: N D^2 constant to 1 % over D 3-7 um at F 0.8', spread < 0.01, 'spread %.2f %%' % (100 * spread))
check('A R13(b) (eq. 82, constant w_c): N D^2 within 1.5 % of -kappa^2 w_c^2 ln(1-F)/6', off < 0.015,
      '%.0f vs %.0f (%.2f %%)' % (np.mean(ND2), closed, 100 * off))
# ---- The exact closures (spec v5 section 13) -------------------------------------------------------------------------
# These closures generalize App. A Result 10, eq. 69 (Gifford eqn 4, App. A eq. 64). The planner assigns
# f(i) = g(delta*/w_c), and the closure chooses the side resistances and the outlet so that the resistor network,
# solved independently (by measure.solve3, the solver that judges drawn masks), reproduces those f and that F.
for name_, sc_ in (('gifford_physics_linear', sc1), ('the same at D 4 um', sc2)):
    Rs_ = np.asarray(sc_.Rs, float); Rg_ = np.asarray(sc_.Rg, float)
    net = M2.solve3(np.asarray(sc_.Rc, float), Rs_, Rs_, Rg_, Rg_, sc_.Rc_out, sc_.Rs_out, sc_.Rs_out)
    n_ = sc_.N - 1
    dev = np.max(np.abs(net['fA'][:n_] / np.asarray(sc_.f[:n_]) - 1))
    check('gap-drop closure (spec v5 sec. 13; gen. A R10 eq. 69): network gives every designed f (%s)' % name_, dev < 1e-8,
          '%.1e' % dev)
    D_ = S.KAPPA_G * sc_.spec.delta_star
    cut = S.KAPPA_G * np.array([S.x_of(f) for f in net['fA'][:n_]]) * np.asarray(sc_.wc[:n_])
    check('  ... and with A R2 and R6, every gap of that network realises D (%s)' % name_, np.max(np.abs(cut / D_ - 1)) < 1e-8)
    check('outlet closure (A R8 last para.; App. B P-C4.3): the network delivers the design F (%s)' % name_,
          abs(net['F'] - sc_.F) < 1e-9, '%.1e' % abs(net['F'] - sc_.F))
    # By App. A Result 11 (eqs. 75 and 76), the march itself satisfies the outlet identity with R_s = R_c/sigma(N) (the
    # PE closure, eq. 69), that is, 1 - P(N) = 2 sigma/(1 + 2 sigma) with sigma = S(N)/P(N), to rounding. This is the
    # implementers' step-7 check.
    PN, SN = float(sc_.P[sc_.N - 1]), float(sc_.S[sc_.N - 1])
    sg = SN / PN
    e76 = abs((1 - PN) - 2 * sg / (1 + 2 * sg))
    check('A R11 (eq. 75): the march meets 1 - P(N) = 2 sigma/(1 + 2 sigma) to rounding (%s)' % name_, e76 < 1e-12, '%.1e' % e76)
    # opencif closes its outlet on the exact gap-drop field rather than at pressure equilibration. Thus eq. 44,
    # evaluated on its own outlet resistances, differs from F by the outlet row's gap drop R_gap Q_gap, which is a few
    # f_gap and would be zero without that drop.
    Fid = 2 * sc_.Rc_out / (sc_.Rs_out + 2 * sc_.Rc_out)
    fm = float(np.max(sc_.f))
    check('exact outlet vs A R8 eq. 44: |F_identity - F| <= 4 f_gap, the outlet row gap drop (%s)' % name_,
          abs(Fid - sc_.F) <= 4 * fm, '%.1f f_gap' % (abs(Fid - sc_.F) / fm))
check('A R10: the design F is the target within one gap fraction', abs(sc1.F - F_) < 2 * np.max(sc1.f))
d3 = copy.deepcopy(base); d3['targets']['f_gap'] = 3e-4
try:
    SP.load(d3); ok3 = False
except SP.SpecError:
    ok3 = True
check('App. B thesis: inputs are D and F; f_gap is internal (a spec giving f_gap is refused)', ok3)
# ---- App. B C7 (P-C7.1): the kappa(u) Hill curve below the knee and Dinh's band --------------------------------------
check('glossary \\hillKmax, \\hillUhalf, \\hillN = schedule.HILL (%s, %s, %s)' % (mac.get('hillKmax'), mac.get('hillUhalf'),
      mac.get('hillN')), S.HILL == (float(mac['hillKmax']), float(mac['hillUhalf']), float(mac['hillN'])))
band, kl = [], None
for wv in (75.0, 52.0):
    dl = wv * math.sqrt(1.04e-4 / 3.0)                       # App. B P-C7.1: delta = w_c sqrt(f/3), f held at 1.04e-4
    k = float(BU.hill_kappa(S.KAPPA_G * dl / wv))
    kl = kl or k
    band += [k * dl, S.KAPPA_G * dl]
check('App. B P-C7.1: Dinh band [\\dinhBandLo, \\dinhBandHi] = [%s, %s] um' % (mac['dinhBandLo'], mac['dinhBandHi']),
      round(min(band), 2) == float(mac['dinhBandLo']) and round(max(band), 2) == float(mac['dinhBandHi']),
      '[%.3f, %.3f]' % (min(band), max(band)))
check('App. B P-C7.1: local kappa \\dinhKappaLocal = %s; Xia-scaled collapse edge \\dinhCollapseLo = %s um'
      % (mac['dinhKappaLocal'], mac['dinhCollapseLo']),
      round(kl, 2) == float(mac['dinhKappaLocal']) and round(kl * 75.0 * math.sqrt(1.04e-4 / 3.0), 2) == float(mac['dinhCollapseLo']))
check('glossary \\kappOverRootThree = kappa_g/sqrt(3) (%s)' % mac['kappOverRootThree'],
      round(S.KAPPA_G / math.sqrt(3), 4) == float(mac['kappOverRootThree']))
# The gates use that curve, so a constant-delta* device at Dinh's width reports the collapse edge as its passage edge.
import types
wcs = np.linspace(75.0, 52.0, 50)
dst = 2.174 / S.KAPPA_G
ns = types.SimpleNamespace(wc=wcs, P=np.linspace(1, 0.5, 50), f=np.array([S.f_of(dst, w) for w in wcs]),
                           Qg=np.full(50, 1e-4), wL=np.full(50, 19.0), delta=np.full(50, dst), Rg=np.ones(50), Rc=np.ones(50))
sp_ = types.SimpleNamespace(Q_in_uL_min=30.0, h=140.0, wc0=75.0, G=19.0, g_p=27.0, Lam_L=46.0, delta_star=dst)
g_ = BU.gates(sp_, ns, 2.174)
pe = g_['cutoff_band_um']['passage_edge_below_knee']
# The cutoff D = 2.174 um (\dcritGGIn, the cubic-pair reading of Dinh's f* at 75 um) and App. B's 2.170 um
# (kappa_g x 0.442) are 0.2% apart.
check('build.gates below the knee: passage edge = kappa_Hill(u) delta* (\\dinhCollapseLo 1.80 um at Dinh V1, +-0.01)',
      pe is not None and abs(pe - float(mac['dinhCollapseLo'])) < 0.01 and
      abs(g_['cutoff_band_um']['band_um'][1] - S.KAPPA_G * dst) < 1e-12, '%.4f um' % (pe or 0))
ns2 = types.SimpleNamespace(**{k: v for k, v in vars(ns).items()})
ns2.wc = np.linspace(400.0, 300.0, 50)
g2 = BU.gates(sp_, ns2, 2.174)
check('build.gates below the foot (u < 0.012): no band, no passage edge', g2['plateau_u_min'] < S.U_FF and
      g2['cutoff_band_um']['band_um'] is None and g2['cutoff_band_um']['passage_edge_below_knee'] is None)
# ---- App. B thesis scope and CIF_v2 M14: the Stokes gates equal CIF_v2's literal formulas on a no-head device --------
import hashlib
check('reference/CIF_v2.py is the frozen v2.1.0 file (md5 0e0961a5..., the project copy of CIF_v2.py)',
      hashlib.md5(open(os.path.join(ROOT, 'reference', 'CIF_v2.py'), 'rb').read()).hexdigest()
      == '0e0961a50441afb94671b4470d41bace' and C.__version__ == '2.1.0')
for q_ in (25.0, 28.5, 30.0):
    dq = copy.deepcopy(base); dq.setdefault('operation', {})['flow_uL_min'] = q_
    spq = SP.to_engine(dq)[0]
    scq = S.plan(spq)
    gq = BU.gates(spq, scq, S.KAPPA_G * spq.delta_star)
    pr = C.Process(h=spq.h * S.UM)
    Qin = q_ * 1e-9 / 60.0
    f1 = S.f_of(spq.delta_star, spq.wc0)
    re_ch = C.reynolds(Qin, spq.wc0 * S.UM, spq.h * S.UM, pr)
    re_gap = C.reynolds(f1 * Qin, spq.G * S.UM, spq.h * S.UM, pr)
    a0, b0 = max(spq.wc0, spq.h) * S.UM, min(spq.wc0, spq.h) * S.UM
    re_d = pr.rho_f * C.C_shear(a0 / b0) * Qin / (a0 * b0 ** 2) * (spq.delta_star * S.UM) ** 2 / pr.mu
    dv = max(abs(gq['Re_max'] / re_ch - 1), abs(gq['Re_gap_max'] / re_gap - 1), abs(gq['Re_delta_max'] / re_d - 1))
    check('CIF_v2 M14 literal (row 1, full inlet flow) = build.gates at %g uL/min: Re %.3f, Re_gap %.3e, Re_delta %.3e'
          % (q_, re_ch, re_gap, re_d), dv < 1e-9 and gq['stokes_ok'] == (re_gap <= 1e-2 and re_d <= 1e-2), '%.1e' % dv)
# ---- CIF_v2 E11 and M16 (march()) against build.gates on matched requests just below and just above each limit -------
# E11 uses a constant w_c of 60 um, h = 20 um, and delta*/h = 0.10 (1 -+ 1e-6). M16 uses a linear W(t) on CIF_v2's
# default process (Lambda = G + g_p = 20.5 + 25.6 = 46.1 um), with steps of 0.45 and 0.47 um per row against
# eps_taper Lambda = 0.461 um. CIF_v2's march() either raises GateError or does not, build.gates' lamina_ok and
# taper_ok must agree with it, and the taper rate must be CIF_v2's number.
UM_ = C.UM


def v2_march(h_um, w_of_t, N, dstar_um, F_t):
    """Run CIF_v2's march() for N rows on the width profile w_of_t at lamina dstar_um and depth h_um, and return its
    gate error (None when it raises none), its rows, and its process."""
    pr = C.Process(h=h_um * UM_)
    req = C.Request(D_target=C.Pairing().kappa * dstar_um * UM_, F_target=F_t, W_shape=lambda t: w_of_t(t) * UM_, proc=pr)
    try:
        rows_, _ = C.march(req, N, dstar_um * UM_)
        return None, rows_, pr
    except C.GateError as e:
        return str(e), None, pr


def ob_gates(h_um, G_um, gp_um, wc0_um, wcs, dstar_um):
    """Evaluate opencif's build.gates on a stand-in schedule with the center widths wcs and a constant lamina dstar_um,
    so that the gates can be compared with CIF_v2's march() on a matched request."""
    wcs = np.asarray(wcs, float)
    n_ = len(wcs)
    ns_ = types.SimpleNamespace(wc=wcs, P=np.linspace(1, 0.9, n_), f=np.array([S.f_of(dstar_um, w) for w in wcs]),
                                Qg=np.full(n_, 1e-4), wL=np.full(n_, G_um), delta=np.full(n_, dstar_um),
                                Rg=np.ones(n_), Rc=np.ones(n_))
    sp2 = types.SimpleNamespace(Q_in_uL_min=1.0, h=h_um, wc0=wc0_um, G=G_um, g_p=gp_um, Lam_L=G_um + gp_um,
                                delta_star=dstar_um)
    return BU.gates(sp2, ns_, S.KAPPA_G * dstar_um)


e11 = []
for fac in (1 - 1e-6, 1 + 1e-6):
    ds = 0.10 * 20.0 * fac
    err, _, _ = v2_march(20.0, lambda t: 60.0, 40, ds, 0.02)
    gg = ob_gates(20.0, 19.0, 30.0, 60.0, [60.0] * 5, ds)
    e11.append(((err is not None and 'eps_lam' in err), (not gg['lamina_ok'])))
check('CIF_v2 E11 (march, eps_lam 0.10) = build.gates lamina_ok at delta*/h = 0.10 -+ 1e-6 (no error, error)',
      e11 == [(False, False), (True, True)], str(e11))
m16 = []
for step in (0.45, 0.47):
    N_ = 40
    wf = lambda t, s_=step, n_=N_: 60.0 - s_ * n_ * t
    err, rows_, pr = v2_march(140.0, wf, N_, 1.0, 0.02)
    lam_um = pr.Lambda / UM_
    n_rows = len(rows_) if rows_ is not None else 1
    gg = ob_gates(140.0, pr.G / UM_, pr.g_p / UM_, wf(0.0), [wf(i / N_) for i in range(1, max(n_rows, 2) + 1)], 1.0)
    v2_rate = float(re.search(r'= ([0-9.eE+-]+) > eps_taper', err).group(1)) if err else None
    m16.append(((err is not None), (not gg['taper_ok']), abs(gg['taper_rate_max'] - step / lam_um) < 1e-12,
                (v2_rate is None or abs(v2_rate - gg['taper_rate_max']) < 1e-3 * v2_rate)))
check('CIF_v2 M16 (march, |w_c(i) - w_c(i-1)|/Lambda, eps_taper 1e-2) = build.gates taper_ok and its rate (steps 0.45, '
      '0.47 um on Lambda 46.1 um)', m16 == [(False, False, True, True), (True, True, True, True)], str(m16))
# The floor w_min = max(F_min, h/AR_max) of opencif (spec.to_engine) is compared with CIF_v2's Process.w_min at
# shallow and deep processes (h = 100, 300, and 80 um).
wm = []
for depth_, floor_ in ((100.0, 19.0), (300.0, 19.0), (80.0, 4.5)):
    dw = copy.deepcopy(base); dw['process']['depth_um'] = depth_; dw['process']['min_feature_um'] = floor_
    dw['process'].pop('aspect_ratio_max', None)                  # the default (the example declares 25)
    o_ = SP.to_engine(dw)[1]
    wm.append(abs(o_['w_min'] - C.Process(h=depth_ * UM_, F_min=floor_ * UM_).w_min / UM_) < 1e-9 and o_['AR_max'] == 10.0)
check('CIF_v2 Process.w_min = max(F_min, h/AR_max) = opencif opts w_min (AR_max default 10) at h 100, 300, 80 um', all(wm))
# CIF_v2's validator G7 floors every w_c at w_min, and wmin_check refuses a planned center width below it (strict mode).
spw, opw = SP.to_engine(copy.deepcopy(base))
opw = dict(opw, mode='strict', w_min=16.0, AR_max=5.0); opw.pop('preflight_refusal', None)
nsw = types.SimpleNamespace(ws=np.array([np.nan, 40.0, 60.0]), wc=np.array([15.0, 15.0, 15.0]))
PL.wmin_check(spw, opw, [], nsw)
check('CIF_v2 validator G7 (every w_c >= w_min): a planned 15 um centre channel under w_min 16 um is refused (strict)',
      'centre width w_c = 15.00 um below' in (opw.get('preflight_refusal') or ''), (opw.get('preflight_refusal') or '')[:60])
# Spec v5 G7 ("widths >= w_min") also applies to the straight side channels. G11 would fold them, whereas opencif
# refuses them in strict mode.
opw2 = dict(opw, mode='strict', w_min=16.0, AR_max=5.0); opw2.pop('preflight_refusal', None)
nsw2 = types.SimpleNamespace(ws=np.array([12.0, 40.0, 60.0, np.nan]), wc=np.array([20.0, 20.0, 20.0, 20.0]), NM=0, N=4)
PL.wmin_check(spw, opw2, [], nsw2)
opw3 = dict(opw, mode='strict', w_min=16.0, AR_max=5.0); opw3.pop('preflight_refusal', None)
nsw3 = types.SimpleNamespace(ws=np.array([12.0, 40.0, 60.0, np.nan]), wc=np.array([20.0, 20.0, 20.0, 20.0]), NM=1, N=4)
PL.wmin_check(spw, opw3, [], nsw3)
check('spec v5 G7 on the straight side channels (G11 would fold them): a planned 12 um straight side channel under w_min '
      '16 um is refused (strict); the same width as a head row (row < NM) is not a straight side channel',
      'straight side channel w_s = 12.00 um (row 1) below' in (opw2.get('preflight_refusal') or '')
      and not opw3.get('preflight_refusal'), (opw2.get('preflight_refusal') or '')[:60])
# The Dean term is checked at the Re = 10 edge, because the first rc3 build compared without a tolerance and dropped
# the term there.
spd = types.SimpleNamespace(Q_in_uL_min=30.0, h=100.0, wc0=100.0, G=19.0, g_p=27.0, Lam_L=46.0, delta_star=1.0, F_min=0.1)
nsd = types.SimpleNamespace(wc=np.full(3, 100.0), P=np.linspace(1, 0.9, 3), f=np.full(3, 1e-4), Qg=np.full(3, 1e-4),
                            wL=np.full(3, 19.0), delta=np.full(3, 1.0), Rg=np.ones(3), Rc=np.ones(3))
gde = BU.gates(spd, nsd, S.KAPPA_G)
check('Dean term given at the edge (Re bound 2 rho Q/(mu h) = 10 to rounding: h 100 um, 30 uL/min, floor 0) and its value',
      gde['dean_leading_term_bound'] is not None and abs(gde['Re_any_channel_bound'] - 10.0) < 1e-9
      and abs(gde['dean_leading_term_bound'] - 0.03058 * (200.0 / 576.0) ** 2) < 1e-12, '%.3e' % (gde['dean_leading_term_bound'] or 0))
# P8 is CIF_v2's embedded physics auditor. The checks below read its contract from CIF_v2's own source text and test
# the record's recomputation of it.
src_p8 = C._PHYSICS_SRC
m_re = re.search(r'Re_max: float = ([0-9.]+)', src_p8)
m_q = re.search(r'Q_in_uLmin: float = ([0-9.]+)', src_p8)
check('P8 (CIF_v2 _PHYSICS_SRC): inlet Re <= Re_max = %s at Q_in_uLmin = %s, and L_e = 0.05 Re D_h <= Lambda/2; build matches'
      % (m_re.group(1) if m_re else '?', m_q.group(1) if m_q else '?'),
      m_re is not None and float(m_re.group(1)) == BU.P8_RE_MAX and m_q is not None and float(m_q.group(1)) == BU.P8_Q_UL_MIN
      and 'Le = 0.05 * Re * Dh' in src_p8 and 'Le <= 0.5 * geom["Lambda"]' in src_p8 and 'ubar = Q_SI / (wc0_SI * h_SI)' in src_p8)
check('CIF_v2 entrance_length (S16-c) itself calls 0.05 Re D_h "simply the wrong formula" (why P8\'s L_e is not gated)',
      'simply the wrong formula' in inspect.getsource(C.entrance_length))
gq8 = BU.gates(spq, scq, S.KAPPA_G * spq.delta_star)                   # the last flow of the M14 loop above (30 uL/min)
pr8 = C.Process(h=spq.h * S.UM)
p8d, p8a = gq8['P8']['at_design_flow'], gq8['P8']['at_auditor_flow']
check('P8 inlet Re in the record = CIF_v2 reynolds(Q, w_c0, h) at the design flow (30) and the auditor\'s (5 uL/min); '
      'the any-channel bound exceeds the centre-channel Re',
      abs(p8d['inlet_Re'] / C.reynolds(30e-9 / 60.0, spq.wc0 * S.UM, spq.h * S.UM, pr8) - 1) < 1e-12
      and abs(p8a['inlet_Re'] / C.reynolds(5e-9 / 60.0, spq.wc0 * S.UM, spq.h * S.UM, pr8) - 1) < 1e-12
      and gq8['Re_any_channel_bound'] > gq8['Re_max'], '%.3f, %.3f' % (p8d['inlet_Re'], p8a['inlet_Re']))
# ---- The numeric macros of the glossary that the code uses -----------------------------------------------------------
check('glossary \\kappVal = kappa_g (%s)' % mac.get('kappVal'), float(mac['kappVal']) == S.KAPPA_G)
check('glossary \\kappCIlo, \\kappCIhi = its interval to 2 d.p. (%s, %s)' % (mac.get('kappCIlo'), mac.get('kappCIhi')),
      round(S.KAPPA_CI[0], 2) == float(mac['kappCIlo']) and round(S.KAPPA_CI[1], 2) == float(mac['kappCIhi']))
check('glossary \\plateauUwc, \\footUwc = the envelope edges the gates use (%s, %s)' % (mac.get('plateauUwc'), mac.get('footUwc')),
      float(mac['plateauUwc']) == S.U_LO and float(mac['footUwc']) == S.U_FF)
cap_ = (re.search(r'\\newcommand\{\\fgapDesignCap\}\{([0-9.]+)\\times10\^\{(-?\d+)\}\}', pre).groups()
        if pre is not None else GLOSSARY['fgapDesignCap'])
check('glossary \\fgapDesignCap = the calibration cap the record reports (schedule.FGAP_CAP)',
      abs(float(cap_[0]) * 10 ** int(cap_[1]) - S.FGAP_CAP) < 1e-12)
import yaml
dv = yaml.safe_load(open(os.path.join(ROOT, 'presets', 'dinh2024_v1_specified.yaml')))
check('glossary \\dcritGGIn, \\FclosedVone = the Dinh V1 preset targets D, F (%s, %s)' % (mac.get('dcritGGIn'), mac.get('FclosedVone')),
      round(dv['targets']['cutoff_um'], 3) == float(mac['dcritGGIn']) and dv['targets']['filtrate_fraction'] == float(mac['FclosedVone']))
from opencif import published as published
check('glossary \\FclosedVthree = published V3 closed-form F (%s)' % mac.get('FclosedVthree'),
      published.Fcl['V3'] == float(mac['FclosedVthree']))
# ---- CIF_v2 correspondence: every supplied constant and tolerance of the design code is consumed or dispositioned ----
# A constant that is neither would be an unported gate that nobody named, as gamma_GS and f_soft were in the review of
# rc1 (F3).
import dataclasses
lim = open(os.path.join(ROOT, 'docs', 'LIMITATIONS.md')).read()
lim_items = {int(m.group(1)): m.group(2) for m in re.finditer(r'(?ms)^(\d+)\. (.*?)(?=^\d+\. |^## |\Z)', lim)}
CONSUMED = 'consumed'                  # a value the code uses (a constant of the law, a band edge)
GATED = 'gated'                        # a gate the code evaluates on every design (and refuses in strict mode)
PARTLY = 'partly gated'                # gated on part of what CIF_v2 gates (gamma_GS only when a head is
                                       # drawn), with the rest reported or not applicable (ALIGNMENT AG7c, AG7f)
REPORTED = 'reported'                  # evaluated on every design and reported, and never a pass condition
DISPOSITION = {
    'Pairing': dict(kernel=(CONSUMED, BU.KERNEL_PAIR['kernel'] == C.Pairing().kernel),
                    kappa=(CONSUMED, S.KAPPA_G == C.Pairing().kappa),
                    kappa_ci=(CONSUMED, tuple(S.KAPPA_CI) == tuple(C.Pairing().kappa_ci)),
                    beta_cal=(REPORTED, tuple(DCR.K.Pairing().beta_cal) == tuple(C.Pairing().beta_cal)),
                    # The transfer bias that beta_cal implies is reported on every build (design_code B7) and never
                    # applied, and on the cubic pair it does not enter N (glossary; ALIGNMENT AG2).
                    u_lo=(GATED, S.U_LO == C.Pairing().u_lo), u_ff=(GATED, S.U_FF == C.Pairing().u_ff),
                    eta=(REPORTED, DCR.ETA == C.Pairing().eta and 'η' in lim_items.get(5, '')),   # eta F, reported (B6)
                    f_soft=(REPORTED, BU.F_SOFT == C.Pairing().f_soft),
                    # gamma_GS is gated when a head is drawn, whereas CIF_v2 gates it always.
                    gamma_GS=(PARTLY, BU.GAMMA_GS == tuple(C.Pairing().gamma_GS))),
    'Tolerances': dict(Re_max=(GATED, BU.RE_MAX == C.Tolerances().Re_max),
                       Re_stokes=(GATED, BU.RE_STOKES == C.Tolerances().Re_stokes),
                       gap_dom=(GATED, BU.M13_MAX == C.Tolerances().gap_dom),
                       eps_lam=(GATED, BU.EPS_LAM == C.Tolerances().eps_lam),
                       eps_taper=(GATED, BU.EPS_TAPER == C.Tolerances().eps_taper),
                       **{k: ('CIF_v2 solver and tripwire tolerance; opencif solves the closures exactly and checks them '
                              'with its own tolerances (this file, Q1/Q2)', True)
                          for k in ('eps_r', 'eps_inv', 'eps_cub', 'eps_id', 'eps_geo', 'eps_net')}),
    'Process': dict(**{k: ('spec input (process / meander / operation / outlet sections)', True)
                       for k in ('h', 'F_min', 'R_min', 'G', 'g_p', 'l_g', 'c_ent', 'GS', 'W', 'Q_in', 'L_out')},
                    mu=('LIMITATIONS 2 (water)', 'water' in lim_items.get(2, '').lower()),
                    rho_f=('LIMITATIONS 2 (water)', 'water' in lim_items.get(2, '').lower()),
                    AR_max=(GATED, BU.AR_MAX_DEFAULT == C.Process().AR_max == SP.DEFAULTS[('process', 'aspect_ratio_max')]
                            and 'AR_max' in lim_items.get(7, '')),   # w_min on G, GS, w_c, and every straight w_s
                                                                     # before drawing and on the drawn channels after
                                                                     # (strict), with the refusal in place of G11's
                                                                     # fold (LIMITATIONS 7)
                    L_eq=('superseded by the exact curved-duct solution for arcs and, for square turns, the computed '
                          'corner term c_corner (ALIGNMENT AG7b)', True),
                    l_stroke=('CIF_v2 fold stroke; replaced by opencif\'s head styles (meander.*)', True),
                    n_legs=('layout.legs', True), leg_length=('chip area: reported, not a pass condition (App. B C8 CQ5)', True),
                    m_budget=('CIF_v2 terminal-layer allowance; opencif closes the outlet exactly', True)),
}
for cls, table in DISPOSITION.items():
    names = {f.name for f in dataclasses.fields(getattr(C, cls))}
    n_ = {lab: sum(1 for how, _ in table.values() if how == lab) for lab in (GATED, PARTLY, REPORTED, CONSUMED)}
    check('CIF_v2 %s: every field accounted for (%d: %d gated, %d partly, %d reported, %d consumed, %d a spec input, '
          'replaced or a named limitation)' % (cls, len(names), n_[GATED], n_[PARTLY], n_[REPORTED], n_[CONSUMED],
                                                len(table) - sum(n_.values())), names == set(table),
          'unmapped %s; stale %s' % (sorted(names - set(table)), sorted(set(table) - names)))
    bad = [k for k, (how, ok) in table.items() if not ok]
    check('CIF_v2 %s: each disposition holds (values equal, or the named LIMITATIONS item names it)' % cls, not bad, str(bad))
# ---- The frozen foundation of the alignment audit: Appendices A and B and the glossary at v114 -----------------------
FROZEN_V114 = {'APPENDIX_A_revised.tex': 'ab28dc46e4e184bc1b0c6d689a1ee129', 'APPENDIX_B_current.tex': '69b565d0cb2ef73a12778cb2df677732',
               'glossary_official.tex': 'e9dc24fdead748e9df9e9eebd608a62e', 'PREAMBLE_FINAL.tex': '3c075af8f63bb1d735f4b43c64f32d14'}
md5 = lambda p: hashlib.md5(open(p, 'rb').read()).hexdigest()
check('reference/glossary_values.json names the frozen v114 preamble (md5 3c075af8...)',
      GLOSSARY['preamble_md5'] == FROZEN_V114['PREAMBLE_FINAL.tex'])
if pre is not None:                     # the dissertation's overlay
    check('reference_glossary/PREAMBLE_FINAL.tex is the frozen v114 preamble (md5 3c075af8...)',
          md5(PRE_PATH) == FROZEN_V114['PREAMBLE_FINAL.tex'])
    m_ = {m.group(1): m.group(2) for m in re.finditer(r'\\newcommand\{\\(\w+)\}\{([0-9.]+)\}', pre)}
    check('reference/glossary_values.json holds the preamble\'s own values (%d macros)' % len(GLOSSARY['macros']),
          all(m_.get(k) == v for k, v in GLOSSARY['macros'].items()) and list(cap_) == GLOSSARY['fgapDesignCap'])
fdir = os.path.join(ROOT, 'foundation')
if os.path.isdir(fdir):                                   # the evaluation archive ships the frozen texts
    got = {f: md5(os.path.join(fdir, f)) for f in FROZEN_V114 if os.path.exists(os.path.join(fdir, f))}
    check('foundation/ holds the frozen v114 A, B, glossary and preamble (md5 as FREEZE_v114)', got == FROZEN_V114,
          str({f: v[:8] for f, v in got.items() if v != FROZEN_V114.get(f)}))
    al = open(os.path.join(ROOT, 'docs', 'ALIGNMENT.md')).read()
    check('docs/ALIGNMENT.md says it was audited against v114', 'frozen v114 sources' in al)
print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
