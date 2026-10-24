"""The code cells that the 1.0.0 tutorial adds to the earlier tutorial: Figure 1, Gifford's Fig. 2 from the planner,
and the solutions of the exercises. make_tutorial_v1.py inserts each one where its section needs it.
"""

FIG1_CELL = r'''
# Figure 1 is drawn by tools/tutorial_figures.py with matplotlib alone, so running it here also shows that the figure
# can be regenerated from the repository.
p_fig = subprocess.run([sys.executable, str(TREE / 'tools' / 'tutorial_figures.py'), str(OUT / 'notebook_figures')],
                       capture_output=True, text=True)
if p_fig.returncode != 0:
    raise RuntimeError('tools/tutorial_figures.py failed:\n' + p_fig.stderr)
preview(OUT / 'notebook_figures' / 'lamina_sketch.png')
'''

FIG2_CELL = r'''
# Gifford et al. (2014), Fig. 2, from opencif's planner. Panels A and B use six constant fractions at a depth of 100 um,
# and panels C and D use f_gap = 1.76e-4 at six depths. Each plan asks for F = 2/3, the paper's endpoint w_s = w_c, and
# for the cutoff that f_gap implies at w_c = 100 um, so that the planner's constant f(i) equals f_gap. The process
# floor and the aspect-ratio limit lie far below the first side width, because the curves are hydraulics alone, and
# the gap's entrance term is zero, because Gifford's eqn (4) contains no gap resistance.
WC_G = 100.0
F_AB = [k * 9.6e-5 for k in (6, 5, 4, 3, 2, 1)]                   # the caption's six values, in equal steps
H_CD = [150.0, 100.0, 50.0, 25.0, 10.0, 5.0]
F_CD = 1.76e-4


def fig2_plan(f, h):
    D = S.KAPPA_G * S.x_of(f) * WC_G                                # the cutoff that f implies (section 3.3, inverted)
    d = opencif.build_spec(D, 2.0 / 3.0, process=dict(depth_um=h, gap_um=20.5, post_length_um=25.6,
                                                      post_width_um=50.6, min_feature_um=0.2,
                                                      aspect_ratio_max=5000.0, gap_entrance_c=0.0),
                           width_in_um=WC_G, meander=None, mode='evaluate', name='fig2')
    S.clear_caches()
    pl = S.plan(SP.to_engine(SP.load(d))[0])
    i = np.arange(1, pl.N + 1, dtype=float)
    return dict(N=pl.N, iT=math.log(3.0) / -math.log(1.0 - 2.0 * f), i=i, sigma=pl.S / pl.P,
                closed=((1.0 - 2.0 * f) ** (-i) - 1.0) / 2.0, ws=pl.ws[:pl.N - 1],
                f_dev=float(np.max(np.abs(pl.f / f - 1.0))))


t0 = time.time()
AB = [fig2_plan(f, 100.0) for f in F_AB]
CD = [fig2_plan(F_CD, h) for h in H_CD]
print('twelve plans in %.0f s' % (time.time() - t0))
print('panels A and B, at a depth of 100 um')
print('    f_gap      N      i_T    largest |f(i)/f_gap - 1|   largest |sigma(i) - closed form|   w_s(N-1) (um)')
for f, g in zip(F_AB, AB):
    print('  %.2e  %5d  %7.1f         %.1e                     %.1e                  %6.2f'
          % (f, g['N'], g['iT'], g['f_dev'], float(np.max(np.abs(g['sigma'] - g['closed']))), g['ws'][-1]))
spread = max(float(np.max(np.abs(g['sigma'] - CD[0]['sigma']))) for g in CD)
print('panels C and D, at f_gap = %.2e, give N = %s at all six depths, and their flow-ratio curves differ by at '
      'most %.1e' % (F_CD, ', '.join(str(n) for n in sorted({g['N'] for g in CD})), spread))
print('the side width at row 1000 is %.1f um at a depth of 150 um and %.1f um at 5 um'
      % (CD[0]['ws'][999], CD[-1]['ws'][999]))
# The rows at which the six curves of the published panel A end, read from the figure (tests/test_design_code.py,
# section 6, gives the method).
DIGITIZED = [951, 1143, 1426, 1898, 2840, 5627]
print('the published curves of panel A end near rows %s, and the plans end at rows %s'
      % (', '.join(map(str, DIGITIZED)), ', '.join(str(g['N']) for g in AB)))

blues = ['#0b2e59', '#163f73', '#1f5aa6', '#2a78d6', '#6f9fdc', '#a9c4ea']
fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.4))
for g, f, c in zip(AB, F_AB, blues):
    ax[0, 0].plot(g['i'][:-1], g['ws'], color=c, lw=1.3, label='%.2e' % f)
    ax[0, 1].plot(g['i'], g['sigma'], color=c, lw=1.3)
for g, h, c in zip(CD, H_CD, blues):
    ax[1, 0].plot(g['i'][:-1], g['ws'], color=c, lw=1.3, label='%g um' % h)
    ax[1, 1].plot(g['i'], g['sigma'], color=c, lw=1.3)
ax[0, 0].set_title('(A) Side-channel width for six values of $f_{gap}$ (depth 100 um)', fontsize=9.5, loc='left')
ax[0, 1].set_title('(B) Ratio of side flow to center flow for the same six values', fontsize=9.5, loc='left')
ax[1, 0].set_title('(C) Side-channel width for six depths ($f_{gap}$ = 1.76 × 10$^{-4}$)', fontsize=9.5, loc='left')
ax[1, 1].set_title('(D) Ratio of side flow to center flow for the six depths', fontsize=9.5, loc='left')
ax[0, 0].legend(title='$f_{gap}$', fontsize=7.5, title_fontsize=8)
ax[1, 0].legend(title='depth', fontsize=7.5, title_fontsize=8)
for a in ax.flat:
    a.set_xlabel('gap row $i$')
    a.grid(alpha=0.3)
for a in ax[:, 0]:
    a.set_ylabel('$w_s(i)$ (um)')
for a in ax[:, 1]:
    a.set_ylabel('$Q_s(i)/Q_c(i)$')
fig.suptitle("Gifford et al. (2014), Fig. 2, reproduced from opencif's plans", fontsize=11)
plt.tight_layout()
show(fig, 'fig_4_5_gifford_fig2')
'''

SOLUTIONS_CELL = r'''
# Exercise 1. The row count scales as 1/D^2 at a fixed F, so the plan at D = 4 um predicts the plan at D = 3.5 um.
sol = OUT / 'solutions'
r4 = opencif.generate(4.0, 0.5, out=sol, check_only='plan', name='ex1_D4', **DC)
r35 = opencif.generate(3.5, 0.5, out=sol, check_only='plan', name='ex1_D3p5', **DC)
r35b = opencif.generate(3.5, 0.6, out=sol, check_only='plan', name='ex1_D3p5_F0p6', **DC)
d4, d35, d35b = (r.record['design'] for r in (r4, r35, r35b))
print('exercise 1: N D^2 = %d at D = 4 um predicts N = %d at D = 3.5 um, and the plan gives N = %d (F planned %.5f, '
      'D planned %.4f um)' % (d4['N'] * 16, round(d4['N'] * 16 / 3.5 ** 2), d35['N'], d35['F_planned'], d35['D_um']))
print('            at F = 0.6 the plan gives N = %d and D planned %.4f um, so N moves and the cutoff does not'
      % (d35b['N'], d35b['D_um']))

# Exercise 2. The plateau requires D/w_c >= 0.035 at the widest center width.
w_max = 3.0 / 0.035
for mode in ('strict', 'evaluate'):
    parts = dict(DC)
    parts['taper'] = ('constant', round(w_max + 2.0, 1))
    r = opencif.generate(3.0, 0.5, out=sol, check_only='plan', mode=mode, name='ex2_' + mode, **parts)
    print('exercise 2: widest center width on the plateau %.1f um; at %.1f um, %s mode gives exit %d (%s)'
          % (w_max, w_max + 2.0, mode, r.exit_code, r.verdict))

# Exercise 3. A sealed start without a head on the 19 um process.
try:
    opencif.generate(4.0, 0.5, meander=None, process='generic19', out=sol, check_only='plan', name='ex3')
    print('exercise 3: no error raised')
except Exception as e:
    print('exercise 3: %s:\n   %s' % (type(e).__name__, e))

# Exercise 4. The exact closure against the PE closure on a tapered device. The center channel narrows within each
# segment, so the PE side resistance is the plan's own segment resistance R_c(i) times P(i)/S(i), and the PE width is
# the one whose resistance over the segment length ell(i) equals it.
def closures(path):
    S.clear_caches()
    p = S.plan(SP.to_engine(SP.load(str(path)))[0])
    seg = np.arange(p.N - 1)
    w_pe = np.array([brentq(lambda w: S.rr(w, p.spec.h) * p.ell[i] * S.UM - p.Rc[i] * p.P[i] / p.S[i], 1e-3, 5e3)
                     for i in seg])
    rs_rel = p.Rs[seg] / (p.Rc[seg] * p.P[seg] / p.S[seg]) - 1
    dq = float(np.median(np.abs(np.diff(p.Qg[:p.N]) / p.Qg[:p.N - 1])))     # row-to-row change of the gap flow
    return p, float(np.max(np.abs(p.ws[seg] / w_pe - 1))), rs_rel, dq


M15 = TREE / 'examples' / 'matrix' / 'm15_linear_shear_nomeander_filter_teardrop.yaml'
p15, dw15, rs15, dq15 = closures(M15)
p_lin, dw_lin, rs_lin, dq_lin = closures(TREE / 'examples' / 'gifford_physics_linear.yaml')
S.clear_caches()
g15 = PL.prepare(SP.load(str(M15)), lambda m: None)[3]
print('exercise 4: m15 has N = %d, w_c from %.0f to %.0f um, f from %.2e to %.2e, and Q_gap/Q_in from %.3e to %.3e'
      % (p15.N, p15.wc[0], p15.wc[p15.N - 1], p15.f[0], p15.f[p15.N - 1], p15.Qg[0], p15.Qg[p15.N - 1]))
print('  largest relative difference of the side widths, exact against PE: %.1e on m15, %.1e on the straight device'
      % (dw15, dw_lin))
print('  R_s(exact) / R_s(PE) - 1 on m15: %.1e (median over the rows; %.1e to %.1e)'
      % (np.median(rs15), rs15.min(), rs15.max()))
print('  row-to-row change of the gap flow (median): %.1e on m15, %.1e on the straight device (2 f = %.1e)'
      % (dq15, dq_lin, 2.0 * p_lin.f[5]))
print('  M13 = %.3f on m15, against its gate of 0.05' % g15['M13_gap_dominance_max'])

# Exercise 5. A full check draws the outline in memory and chooses the legs; the realized cutoff needs a build.
r5 = opencif.generate(5.0, 0.7, chip=(25, 75), layout='paperclip', outlet='ring', out=sol, check_only=True,
                      name='ex5')
lay5 = r5.record.get('layout') or {}
print('exercise 5: exit %d (%s); legs %s' % (r5.exit_code, r5.verdict, lay5.get('n_legs', lay5.get('legs'))))

# Exercise 6. The provenance tags of the Xia preset and the reasons for its exit code.
head = [ln for ln in open(TREE / 'presets' / 'xia2016_rs7_aligned.yaml') if ln.startswith('#')]
print('exercise 6: the preset header has %d comment lines; its first lines are' % len(head))
print(''.join(head[:12]).rstrip())
print('the released build of this preset: exit %d, physics audit %s'
      % (FROZEN['xia2016_rs7_aligned']['exit_code'], FROZEN['xia2016_rs7_aligned']['physics_audit']))
'''
