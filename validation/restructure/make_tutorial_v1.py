"""Write notebooks/opencif_tutorial.ipynb (unexecuted) -- the undergraduate tutorial: from a cutoff and a filtrate
fraction to an audited photomask, with opencif 1.0.0 from the repository it sits in.
usage: python make_tutorial_v1.py OUT.ipynb"""
import sys
import nbformat as nbf

cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s.strip('\n')))
code = lambda s: cells.append(nbf.v4.new_code_cell(s.strip('\n')))

md(r"""
# From a cutoff and a filtrate fraction to a photomask
### Controlled Incremental Filtration (CIF) with opencif 1.0.0: a tutorial

A CIF device sorts particles by size. A suspension flows down a **center channel**; at every **gap row** a thin layer
of fluid next to each wall is skimmed through the gaps into two **side channels**. A particle too large to fit in
that skimmed layer stays in the center channel, and the center outlet collects the larger particles (the
*retentate*); the side outlets collect the fluid and the smaller particles (the *filtrate*). Gifford *et al.*
(*Lab Chip* 2014) introduced the design; Xia *et al.* (*Sci. Rep.* 2016) and Dinh *et al.* (*Lab Chip* 2024) built
devices with it.

A CIF device is designed from **two numbers**:

* the **cutoff** $D$ (the glossary's $d_{\mathrm{crit}}$): particles larger than $D$ are retained;
* the **filtrate fraction** $F$: the share of the inlet flow that leaves through the side channels.

The two are *decoupled*: the cutoff is set by how much each gap skims, the filtrate fraction by how many gap rows
there are (Appendix A, Results 8 and 13; Appendix B, claims C4 and C8). This notebook walks from those two numbers
to a photomask (a DXF file a mask shop can print), checks every step, and then audits the photomask itself, read back
from its file.

**What runs here.** opencif 1.0.0, from the repository this notebook sits in (`notebooks/` of the repository; the
code stamp is printed and checked in section 1). Every equation below is cited to the frozen foundation of the
dissertation: **Appendix A** ("Step-by-Step Derivations Behind Controlled Incremental Filtration": Results 1-15),
**Appendix B** (the argumentation: claims C1-C9) and the **glossary**; their texts are in `foundation/`, and the test
suite `tests/test_appendix_trace.py` pins the code to them. Nothing in this notebook adds physics: it calls the
package's functions and prints what they compute; where it computes an equation itself (the march, the row count, the
PE closure) it does so to compare with the package.

**Contents**
1. Setup: load and verify the repository
2. The device and its vocabulary
3. The cutoff axis: from $D$ to the fraction each gap skims
4. The yield axis: from $F$ to the number of rows, and from flows to widths
5. The whole path in 300 lines: a straight ladder by hand (`opencif.teaching`)
6. One call: from $(D, F)$ to a photomask, and what the record checks
7. The physics audit: the photomask read back from its file
8. The validity envelope: what the generator refuses, and why
9. Composing new devices from parts
10. Reproducing the published devices
11. opencif and the design code `CIF_v2.py`
12. Micrograph traces (the companion notebook)
13. Exercises
14. The reproducibility record

**Running it.** Python 3.11 with the tested packages (`requirements-tested.txt`: numpy 2.4.4, scipy 1.17.1, shapely
2.1.2, matplotlib 3.10.9, pyyaml 6.0.3, ezdxf 1.4.4, pillow 12.2.0) and Jupyter, started in the repository's
`notebooks/` folder. Sections 1-8 take about 10-20 minutes; section 9 builds three new devices (about 20 minutes, up
to 50 while section 10's builds run beside them); section 10 builds the four published devices (about 1-1.5 hours on
two cores; set the environment variable `OPENCIF_TUTORIAL_QUICK=1` before starting Jupyter to plan them only, in 3-5
minutes). Outputs go to `notebooks/tutorial_out/`.
""")

md(r"""
## 1. Setup: load and verify the repository

The cell below finds the repository this notebook sits in, imports its `opencif` package, and checks that it is the
released code: every record opencif writes carries a **code stamp**, a sha256 over every Python file of the package
and its data, and `validation/RELEASE.json` holds the stamp of the release and the evidence it was validated with. It
also checks that the design code `reference/CIF_v2.py` is the frozen file Appendix A names. If any file of the package
had been edited, the stamp would differ and the cell would say so.
""")
code(r"""
import os, sys, json, math, time, hashlib, subprocess, shutil, glob
from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from IPython.display import display, Markdown, Image

HERE = Path.cwd()
TREE = next(p for p in [HERE, *HERE.parents] if (p / 'opencif' / '__init__.py').exists())   # the repository
sys.path.insert(0, str(TREE))
import opencif
from opencif import pipeline as PL
RELEASE = json.load(open(TREE / 'validation' / 'RELEASE.json'))      # this release: its code stamp and evidence
OUT = HERE / 'tutorial_out'
OUT.mkdir(exist_ok=True)
QUICK = os.environ.get('OPENCIF_TUTORIAL_QUICK') == '1'
plt.rcParams.update({'figure.dpi': 110, 'savefig.dpi': 110, 'font.size': 9})

env = PL.ENV
md5 = hashlib.md5(open(TREE / 'reference' / 'CIF_v2.py', 'rb').read()).hexdigest()
print('opencif', opencif.__version__, '| code stamp', PL.CODE['sha256'], '| tests', env['tests_sha256'])
print('the released code (validation/RELEASE.json):', RELEASE['code_stamp'],
      '-> this package IS the released code' if PL.CODE['sha256'] == RELEASE['code_stamp'] else
      '-> THIS PACKAGE DIFFERS FROM THE RELEASED CODE: its records will say so')
print('the design code reference/CIF_v2.py: md5 %s (%s)' % (md5, 'the frozen v2.1.0' if md5 == '0e0961a50441afb94671b4470d41bace'
                                                            else 'NOT the frozen v2.1.0'))
print('Python', env['python'], '|', ', '.join('%s %s' % kv for kv in sorted(env['installed'].items())))
print('installed versions equal the tested ones (requirements-tested.txt):', env['matches_lock'])
print('quick mode (plans only for the published devices):', QUICK)
""")
code(r"""
def preview(path, width=1100):
    '''show a PNG scaled to `width` pixels (the file on disk keeps its full resolution)'''
    from PIL import Image as PILImage
    import io
    im = PILImage.open(path)
    if im.width > width:
        im = im.resize((width, round(im.height * width / im.width)), PILImage.LANCZOS)
    buf = io.BytesIO()
    im.convert('RGB').save(buf, format='PNG', optimize=True)
    display(Image(data=buf.getvalue()))
    print('(preview; full resolution: %s)' % os.path.relpath(path, HERE))


def show(fig, name):
    '''display a figure inline as a PNG, and keep that same PNG in tutorial_out/notebook_figures/<name>.png (opencif's
    renderers select matplotlib's Agg backend, under which a Figure object would display only as text)'''
    import io
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    data = buf.getvalue()
    (OUT / 'notebook_figures').mkdir(parents=True, exist_ok=True)
    (OUT / 'notebook_figures' / (name + '.png')).write_bytes(data)
    display(Image(data=data))


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


import threading
REUSED = {}          # reproduction builds found already done by an earlier run (section 14 records them)
REPRO = ['dinh2024_v1_figure_aligned', 'gifford2014_fig4_aligned', 'gifford2014_fig6_aligned', 'xia2016_rs7_aligned']


def build_presets(names, workers=2, log_path=None):
    '''build presets in separate processes, `workers` at a time (python -m opencif build presets/<name>.yaml, from the
    repository's folder); a build already finished here whose photomask's sha256 is the released one is not repeated'''
    say = open(log_path, 'a') if log_path else sys.stdout
    stamp = lambda: time.strftime('%H:%M:%S')
    todo = []
    for n in names:
        dxf = OUT / 'repro' / n / (n + '.dxf')
        if dxf.exists() and (OUT / 'repro' / n / (n + '.json')).exists() and sha256(dxf) == FROZEN[n]['dxf_sha256']:
            rj_ = json.load(open(OUT / 'repro' / n / (n + '.json')))
            REUSED[n] = dict(code_stamp=(rj_.get('code') or {}).get('sha256'), written=time.ctime(os.path.getmtime(dxf)))
            print('%s %-28s already built here (sha256 verified)' % (stamp(), n), file=say, flush=True)
        else:
            todo.append(n)
    running, t0 = {}, time.time()
    while todo or running:
        while todo and len(running) < workers:
            n = todo.pop(0)
            (OUT / 'repro').mkdir(exist_ok=True)
            log = open(OUT / 'repro' / (n + '.log'), 'w')
            running[n] = (subprocess.Popen([sys.executable, '-m', 'opencif', 'build',
                                            str(TREE / 'presets' / (n + '.yaml')), '--out', str(OUT / 'repro' / n)],
                                           cwd=str(TREE), stdout=log, stderr=subprocess.STDOUT), log, time.time())
            print('%s %-28s started' % (stamp(), n), file=say, flush=True)
        for n, (p, log, ts) in list(running.items()):
            if p.poll() is not None:
                log.close()
                print('%s %-28s finished, exit %d, %.0f min' % (stamp(), n, p.returncode, (time.time() - ts) / 60),
                      file=say, flush=True)
                del running[n]
        time.sleep(10)
    print('%s all reproduction builds done: %.0f min' % (stamp(), (time.time() - t0) / 60), file=say, flush=True)


FROZEN = RELEASE['photomasks']     # the release's own rebuild of every shipped spec: photomask sha256, exit code, N, F, audit
print(len(FROZEN), 'released builds on record, e.g. gifford2014_fig4_aligned:',
      {k: FROZEN['gifford2014_fig4_aligned'][k] for k in ('exit_code', 'N', 'F_design', 'dxf_sha256', 'physics_audit')})
""")

md(r"""
The four published devices take about an hour to build (section 10). The cell below starts their builds now, in
separate processes two at a time, so they run while you read on; section 10 collects them. (Quick mode does not.)
""")
code(r"""
(OUT / 'repro').mkdir(exist_ok=True)
REPRO_LOG = OUT / 'repro' / 'builds.log'
repro_thread = None
if not QUICK:
    repro_thread = threading.Thread(target=build_presets, args=(REPRO, 2, REPRO_LOG), daemon=True)
    repro_thread.start()
    print('reproduction builds started in the background; progress: %s' % os.path.relpath(REPRO_LOG, HERE))
""")

md(r"""
## 2. The device and its vocabulary

The symbols below are the glossary's (Table A.1, "Geometric Variables"; the "Master Cheat-Sheet"), except where a
row says it is opencif's usage; the code's names are given where they differ. One naming trap is worth learning now: the glossary calls $g_p$, the post's extent
*along* the flow, the "post width", and $\ell_g$, its extent *across* the flow, the "link length"; the YAML keys are
`post_length_um` ($g_p$) and `post_width_um` ($\ell_g$) (`docs/ALIGNMENT.md`, AG5).

| symbol | meaning (glossary) | in the code |
|---|---|---|
| $w_c$ | center-channel width | `centre_channel.width_in_um` (and `width_min_um` for a taper) |
| $w_s(i)$ | side-channel width of segment $i$ (after gap row $i$) | planned: `plan.ws`; drawn: measured on the mask |
| $h$ | channel depth (every channel) | `process.depth_um` |
| $G$ | link (gap) width, along the flow | `process.gap_um` |
| $g_p$ | post width, along the flow | `process.post_length_um` |
| $\ell_g$ | link length: the post's extent across the flow | `process.post_width_um` |
| $\Lambda \ge G + g_p$ | row pitch (equal on these straight rows; for the glossary a design constraint, not a definition) | |
| $N$ | number of gap rows | `design.N` in the record |
| $f(i)$, $f_{\mathrm{gap}}$ | per-side, per-gap extraction fraction | `plan.f` |
| $\delta$, $x=\delta/w_c$ | diverted-lamina width; reduced skim width | `plan.delta` |
| $\delta^*$, $f^*$ | opencif's usage: the design lamina $D/\kappa_g$, and the per-gap fraction that gives the cutoff at gap 1, $g(\delta^*/w_c)$ (the glossary's $f^*_{\mathrm{gap}}$ is Gifford's empirical threshold, read from data) | `spec.delta_star` |
| $D$ ($d_{\mathrm{crit}}$) | cutoff (critical) diameter | `targets.cutoff_um` |
| $\tilde d = d_{\mathrm{crit}}/w_c$ | reduced size ratio (the records call it $u$) | `validation.gates.plateau_u_min` |
| $\kappa_{\mathrm{app}} \equiv \kappa_g$ | apparent cutoff prefactor, cubic branch: 4.9138 | `opencif.physics.kernels.Pairing().kappa` |
| $F$ | filtrate fraction | `targets.filtrate_fraction` |
| $P(i)$, $S(i)$ | center and one side's flow after row $i$, over the inlet center flow (App. A, Result 10) | `plan.P`, `plan.S` |
| $\mathrm{Re}_{\mathrm{gap}}$, $\mathrm{Re}_\delta$ | opencif's usage (after the design code's M14): Reynolds numbers of the flow through a gap and of the skimmed lamina (the Stokes gates, section 8) | `validation.gates` |
| DRC | opencif's usage: design-rule check, every drawn feature at or above the process floor | `validation.drc` |

The sketch below is drawn here (not to scale; widths exaggerated): five gap rows of a straight device.
""")
code(r"""
from matplotlib.patches import FancyBboxPatch, Rectangle
fig, ax = plt.subplots(figsize=(11, 3.9))
G, gp, wc, lg = 1.0, 1.4, 3.0, 0.9                  # gap, post (along the flow), center width, post (across the flow)
pitch, x0, nrow = G + gp, 1.0, 5                     # gap row i spans [x0 + (i-1) pitch, x0 + (i-1) pitch + G]
L = x0 + nrow * pitch + 0.6                          # drawn length
ws = [0.0, 0.25, 0.45, 0.62, 0.78, 0.92]             # side width after row i (exaggerated; 0 before row 1: a sealed start)
blue, grey = '#cfe3f7', '#8a8a8a'
ax.add_patch(Rectangle((0, -wc / 2), L, wc, color=blue))                       # center channel
for side in (+1, -1):
    y_wall = side * wc / 2                                                       # the center channel's edge
    y_lo = y_wall if side > 0 else y_wall - lg                                   # the post row's band
    ax.add_patch(Rectangle((x0, y_lo), L - x0, lg, color=blue))                  # fluid in the gaps
    ax.add_patch(Rectangle((0, y_lo), x0, lg, color=grey))                       # solid wall upstream of gap 1
    for i in range(nrow):                                                        # a post after each gap
        xp = x0 + i * pitch + G
        ax.add_patch(FancyBboxPatch((xp, y_lo), gp if i < nrow - 1 else L - xp, lg,
                                    boxstyle='round,pad=0,rounding_size=0.2', color=grey))
    for i in range(1, nrow + 1):                                                 # side segment i, after row i
        xa = x0 + (i - 1) * pitch
        xb = x0 + i * pitch if i < nrow else L
        y_s = y_wall + side * lg
        ax.add_patch(Rectangle((xa, y_s if side > 0 else y_s - ws[i]), xb - xa, ws[i], color=blue))
for i in range(1, nrow + 1):
    ax.text(x0 + (i - 1) * pitch + G / 2, -wc / 2 - lg / 2, str(i), ha='center', va='center', fontsize=9)
ax.annotate('', xy=(L - 0.8, 0), xytext=(0.2, 0), arrowprops=dict(arrowstyle='->', color='#2a78d6', lw=1.5))
ax.text(0.3, 0.3, 'center flow $Q_c(i)$', color='#2a78d6', fontsize=9)
xg3 = x0 + 2 * pitch + G / 2
ax.annotate('', xy=(xg3, wc / 2 + lg + 0.45), xytext=(xg3, wc / 2 - 0.35),
            arrowprops=dict(arrowstyle='->', color='#eb6834', lw=1.4))
ax.text(xg3 + 0.15, wc / 2 - 0.75, '$f(3)\\,Q_c(2)$ skimmed at gap row 3', color='#eb6834', fontsize=8)
xm = x0 + 3.5 * pitch
ax.annotate('', xy=(xm, wc / 2), xytext=(xm, -wc / 2), arrowprops=dict(arrowstyle='<->', color='k', lw=0.7))
ax.text(xm + 0.1, -0.9, '$w_c$', fontsize=10)
xs = x0 + 3.5 * pitch
ax.annotate('', xy=(xs, -wc / 2 - lg - ws[4]), xytext=(xs, -wc / 2 - lg), arrowprops=dict(arrowstyle='<->', color='k', lw=0.7))
ax.text(xs + 0.12, -wc / 2 - lg - 0.75, '$w_s(4)$: widens row by row', fontsize=8)
ax.annotate('', xy=(x0, wc / 2 + lg + 1.05), xytext=(x0 + G, wc / 2 + lg + 1.05), arrowprops=dict(arrowstyle='<->', color='k', lw=0.7))
ax.text(x0 + G / 2, wc / 2 + lg + 1.25, 'gap $G$', fontsize=8, ha='center')
ax.annotate('', xy=(x0 + G, wc / 2 + lg + 1.05), xytext=(x0 + pitch, wc / 2 + lg + 1.05), arrowprops=dict(arrowstyle='<->', color='k', lw=0.7))
ax.text(x0 + G + gp / 2, wc / 2 + lg + 1.25, 'post $g_p$', fontsize=8, ha='center')
ax.text(-0.15, wc / 2 + lg / 2, '$\\ell_g$', fontsize=9, ha='right', va='center')
ax.text(L + 0.3, wc / 2 + lg + 0.5, 'side channel: filtrate', fontsize=9, va='center')
ax.text(L + 0.3, 0, 'center channel: retentate', fontsize=9, va='center')
ax.text(L + 0.3, -wc / 2 - lg - 0.5, 'side channel: filtrate', fontsize=9, va='center')
ax.set_xlim(-0.8, L + 4.2); ax.set_ylim(-wc / 2 - lg - 1.6, wc / 2 + lg + 1.7); ax.set_aspect('equal'); ax.axis('off')
ax.set_title('A CIF ladder: gap rows 1-5, numbered in the lower post row (sketch, not to scale)', fontsize=10)
show(fig, 'fig_2_ladder_sketch')
""")

md(r"""
## 3. The cutoff axis: from $D$ to the fraction each gap skims

### 3.1 The skim cubic (Appendix A, Result 1)

In laminar flow between two walls the velocity profile is a parabola. The fraction of the center flow carried by a
layer of width $\delta$ next to one wall, with $x=\delta/w_c$, is

$$ f = 3x^2 - 2x^3 \qquad\text{(App. A, eq. 10)}. $$

The design code `CIF_v2.py` implements it as `g_cubic`; opencif's kernels (`opencif/physics/kernels.py`) restate the
design code's formulas, and `tests/test_kernels.py` proves that they return the same floating-point numbers bit for
bit (the cell below repeats that check on a grid). Appendix A defines the cubic for $0 \le x \le 1$ (eq. 5); opencif's
`f_of` accepts only $x < 1/2$, because beyond it the laminae of the two walls would meet (App. A, "For implementers",
step 3).
""")
code(r"""
from opencif.physics import schedule as S    # opencif's row schedule
from opencif.physics import kernels as K     # Appendix A's kernels
sys.path.insert(0, str(TREE / 'reference'))
import CIF_v2 as C                           # the design code, byte-frozen (the package never imports it)
import inspect
print(inspect.getsource(K.g_cubic))          # the kernel, with its derivation cited
x = np.linspace(0, 0.5, 201)
f = np.array([K.g_cubic(v) for v in x])
print('largest |g_cubic(x) - (3x^2 - 2x^3)| on [0, 1/2]:', np.max(np.abs(f - (3 * x**2 - 2 * x**3))))
print('opencif\'s g_cubic = the design code\'s, bit for bit, at all 201 points:',
      np.array(f).tobytes() == np.array([C.g_cubic(v) for v in x]).tobytes())
fig, ax = plt.subplots(figsize=(5.2, 3.2))
ax.plot(x, f, color='#2a78d6')
ax.set_xlabel('$x = \\delta / w_c$'); ax.set_ylabel('$f$: fraction of the center flow skimmed (one wall)')
ax.set_title('Result 1: $f = 3x^2 - 2x^3$'); ax.grid(alpha=0.3)
show(fig, 'fig_3_1_skim_cubic')
""")

md(r"""
### 3.2 Inverting it (Appendix A, Result 2)

A gap skims a *small* fraction, $f\sim10^{-4}$. Then $x$ is small and the cubic term is negligible:
$\delta \approx w_c\sqrt{f/3}$ (the leading law, App. A eq. 18; App. B claim C2). The exact root of the cubic (by
Viète's formula, App. A eq. 20) is what opencif uses (`kernels.x_of_f`, found by Brent's method). How large is the
difference? Appendix B's calibration cap for the cutoff law is $f \le 5.76\times10^{-4}$ (glossary
`\fgapDesignCap`); the alignment matrix records "exact root vs leading order $\le 0.47\,\%$" up to it (AA6).
""")
code(r"""
for fv in (1e-5, 1e-4, 5.76e-4):
    xe, xl = S.x_of(fv), math.sqrt(fv / 3)
    print('f = %.2e: exact x = %.6f, leading sqrt(f/3) = %.6f, leading/exact - 1 = %+.3f %%'
          % (fv, xe, xl, 100 * (xl / xe - 1)))
""")

md(r"""
### 3.3 The cutoff law, and what $D$ fixes (App. A Results 5-7 and 13a; App. B claim C3)

The cutoff is proportional to the lamina: a particle whose radius exceeds the lamina's width cannot enter it
(the steric argument, App. A Result 5, $\kappa_{\mathrm{steric}}=2$), and the fitted operating law is

$$ d_{\mathrm{crit}} = \kappa_{\mathrm{app}}\,\delta = \kappa_{\mathrm{app}}\, w_c\sqrt{f/3} . $$

Appendix B (claim C3b) fits the plateau prefactor $\kappa_{\mathrm{app}} = 4.91$ [4.75, 5.08] on Gifford's
calibration data; the glossary's cubic-branch value is $\kappa_g = 4.9138$. Measured cutoffs are about 2.5 times
the steric value ($\kappa_{\mathrm{app}}/\kappa_{\mathrm{steric}} = 2.4569$; App. A Result 5, "From the steric value to
the measured prefactor": the finite depth, a per-gap over-extraction and finite-particle effects): the steric picture
fixes the scaling, the calibration fixes the prefactor, and opencif uses the measured $\kappa_g$. The *pairing rule* (App. A, Result 9;
glossary, "Cutoff-Law Prefactor Decomposition") says a prefactor goes with the kernel it was fitted with: opencif uses
the cubic kernel with $\kappa_g$ at every call site (AA13).

So a target cutoff fixes the lamina once, $\delta^* = D/\kappa_g$, and at every row the fraction is
$f(i) = g\!\left(\delta^*/w_c(i)\right)$: the **invariant-lamina schedule** (App. A Result 13a; the design code's S1).
For a constant center width every row skims the same fraction; a narrower center channel must skim *more* to keep the
same lamina, $f \propto w_c^{-2}$ (Result 7).
""")
code(r"""
print('kappa_g = %.4f (95%% interval %.2f-%.2f: App. B C3b prints [4.75, 5.08])' % (S.KAPPA_G, *S.KAPPA_CI))
D = 4.0
dstar = D / S.KAPPA_G
print('D = %.1f um  ->  lamina delta* = D / kappa_g = %.4f um' % (D, dstar))
for wc in (60, 80, 100):
    print('  w_c = %3d um: f = g(delta*/w_c) = %.3e,  reduced size ratio D/w_c = %.4f' % (wc, S.f_of(dstar, wc), D / wc))

wcs = np.linspace(30, 200, 300)
fig, ax = plt.subplots(figsize=(6.4, 3.6))
for Dv, col in ((2.0, '#9ab6d8'), (4.0, '#2a78d6'), (6.0, '#163f73')):
    ok = Dv / wcs >= 0.035                           # on the kappa plateau
    fv = np.array([S.f_of(Dv / S.KAPPA_G, w) for w in wcs])
    ax.plot(wcs[ok], fv[ok], color=col, label='D = %g um' % Dv)
    ax.plot(wcs[~ok], fv[~ok], color=col, ls=':', lw=1)
ax.axhline(5.76e-4, color='#eb6834', lw=1)
ax.text(150, 6.2e-4, 'calibration cap 5.76e-4', color='#eb6834', fontsize=8)
ax.set_yscale('log'); ax.set_xlabel('center width $w_c$ (um)'); ax.set_ylabel('$f$ per gap')
ax.set_title('What a cutoff asks of each gap (dotted: off the plateau, D/w_c < 0.035)')
ax.legend(fontsize=8); ax.grid(alpha=0.3, which='both')
show(fig, 'fig_3_3_f_per_gap')
""")

md(r"""
### 3.4 Where the law holds: the plateau and the calibration cap

$\kappa_{\mathrm{app}}$ is a constant only on the **plateau**, $\tilde d = d_{\mathrm{crit}}/w_c \gtrsim 0.035$
(glossary `\plateauUwc`; App. B, Notation and C3c). Below the plateau edge the prefactor depends on $\tilde d$ and
the cutoff law becomes implicit; it falls through a knee near $\tilde d \approx 0.03$ to a sub-steric *foot* (a prefactor below the steric value 2),
$\kappa_{\mathrm{app}} \approx 1.34$ for $\tilde d < 0.012$ (the mean of three 1 um-bead points, outside the regime the
dissertation's claims cover). Below the plateau edge opencif reports a cutoff *band* instead of a single cutoff
(App. B claim C7: a Hill curve, a smooth S-shaped fit of $\kappa_{\mathrm{app}}$ against $\tilde d$ through the
knee, P-C7.1, anchored on the plateau as C7's CQ3 says; `docs/ALIGNMENT.md` AB5), and none below the foot. So for a given $D$ the
center channel may be at most $D/0.035$ wide. The default of `opencif.generate` is $w_c = D/0.05$, inside the
plateau. And $\kappa_g$ was fitted on gaps skimming up to $5.76\times10^{-4}$; a design above that cap is reported,
not refused (AG12). Section 6 shows what the generator does when a design leaves the plateau.
""")

md(r"""
## 4. The yield axis: from $F$ to the number of rows, and from flows to widths

### 4.1 The march (Appendix A, Result 10)

Number the gap rows $i=1\dots N$. Divide every flow by the inlet center flow: $P(i)$ is the center flow after row
$i$, $S(i)$ one side channel's flow. At row $i$ each side gap takes $f(i)P(i-1)$, so

$$ P(i) = \big(1-2f(i)\big)P(i-1), \qquad S(i) = S(i-1) + f(i)P(i-1), \qquad P(0)=1,\ S(0)=\sigma_0 . $$

With a sealed inlet (all fluid enters the center channel) $\sigma_0=0$. Mass is conserved:
$P(i)+2S(i) = 1+2\sigma_0$ (the **mass identity**). The filtrate fraction is what has left the center channel:

$$ F = 1-P(N) = 1-\prod_{i=1}^{N}\big(1-2f(i)\big) \qquad\text{(App. A eq. 67; the glossary's } F_{\mathrm{recursion}}\text{)}. $$

These five lines are the whole march. Appendix A's implementers' example: $f=10^{-4}$ at every row, 5000 rows, gives
$P = 0.367843$, $\sigma = S/P = 0.859277$, $F = 0.632157$ (AA14).
""")
code(r"""
def march(f_rows, sigma0=0.0):
    P, Sd = 1.0, sigma0
    Ps, Ss = [], []
    for fi in f_rows:
        P, Sd = (1 - 2 * fi) * P, Sd + fi * P        # (the right side uses the old P)
        Ps.append(P); Ss.append(Sd)
    return np.array(Ps), np.array(Ss)

P, Sd = march([1e-4] * 5000)
print('after 5000 rows: P = %.6f, sigma = S/P = %.6f, F = 1 - P = %.6f' % (P[-1], Sd[-1] / P[-1], 1 - P[-1]))
print('largest violation of the mass identity P + 2S = 1:', np.max(np.abs(P + 2 * Sd - 1)))
""")

md(r"""
### 4.2 How many rows? (App. A, Result 10, "Row count")

For a constant fraction, $F = 1-(1-2f)^N$, so $N = \ln(1-F)/\ln(1-2f)$, rounded up. Appendix A's example:
$F = 0.5$ at $f = 1.04\times10^{-4}$ needs 3332.1, so 3333 rows.
""")
code(r"""
Nc = math.log(1 - 0.5) / math.log(1 - 2 * 1.04e-4)
print('N = ln(1 - F) / ln(1 - 2f) = %.1f  ->  %d rows' % (Nc, math.ceil(Nc)))
""")

md(r"""
### 4.3 The design map: $D$ fixes the lamina, $F$ fixes $N$ (App. A Result 13; App. B claim C8)

Because the cutoff fixes $\delta^*$ (and so every $f(i)$), the filtrate target can only be met by the number of
rows. The Jacobian of the map $(\delta^*, N)\mapsto(d_{\mathrm{crit}}, F)$ is lower-triangular: on the plateau, changing $N$ does not move the cutoff (Result 13c). And since each row skims a fraction proportional to $\delta^{*2}$, $N D^2$ is
nearly constant at fixed $F$ and channel shape: halving the cutoff quadruples the rows (Result 13b).

The cell below asks the generator to **plan** nine devices (the parts of `examples/decoupling/`: a straight leg,
$w_c$ = 80 um, a packed meander head, 100 um deep, 19 um gaps) for three cutoffs and three filtrate fractions.
`check_only='plan'` plans and evaluates the gates without drawing: seconds per design.
""")
code(r"""
DC = dict(process=dict(depth_um=100, gap_um=19, post_length_um=27, post_width_um=50, min_feature_um=19),
          layout='linear', taper=('constant', 80),
          meander=dict(style='pack', channel_um=28, wall_um=27, head_width_um=4000, seal='require'),
          inlet=dict(style='plain', port_radius_um=750),
          outlet=dict(style='separate', port_radius_um=750, centre_port_radius_um=750))
rows = []
for Dv in (3.0, 4.0, 5.0):
    for Fv in (0.3, 0.5, 0.7):
        r = opencif.generate(Dv, Fv, out=OUT / 'plans', check_only='plan', name='dc_D%g_F%g' % (Dv, Fv), **DC)
        d = r.record['design']
        rows.append((Dv, Fv, d['N'], d['F_planned'], d['D_um'], r.exit_code))
print(' D (um)  F target   N     F planned   D planned  exit')
for Dv, Fv, N, Fp, Dp, ec in rows:
    print('  %.1f     %.1f    %5d    %.5f     %.4f     %d' % (Dv, Fv, N, Fp, Dp, ec))
for Fv in (0.3, 0.5, 0.7):
    nd = [N * Dv**2 for Dv, F_, N, *_ in rows if F_ == Fv]
    print('F = %.1f: N D^2 = %s  (spread %.2f %%)' % (Fv, ', '.join('%.0f' % v for v in nd), 100 * (max(nd) / min(nd) - 1)))
fig, ax = plt.subplots(figsize=(5.4, 3.4))
for Fv, col in ((0.3, '#9ab6d8'), (0.5, '#2a78d6'), (0.7, '#163f73')):
    Ds = [Dv for Dv, F_, *_ in rows if F_ == Fv]; Ns = [N for Dv, F_, N, *_ in rows if F_ == Fv]
    ax.loglog(Ds, Ns, 'o-', color=col, label='F = %.1f' % Fv)
ax.loglog([3, 5], [620, 620 * (3 / 5) ** 2], 'k--', lw=0.8, label='slope -2: $N \\propto D^{-2}$')
ax.set_xticks([3, 4, 5]); ax.set_xticklabels(['3', '4', '5']); ax.minorticks_off()
ax.set_yticks([250, 500, 1000, 2000, 4000]); ax.set_yticklabels(['250', '500', '1000', '2000', '4000'])
ax.set_xlabel('cutoff D (um)'); ax.set_ylabel('rows N'); ax.grid(alpha=0.3, which='both')
ax.legend(fontsize=8, loc='center left', bbox_to_anchor=(1.02, 0.5))
ax.set_title('Decoupling: D sets the lamina, F sets N')
show(fig, 'fig_4_3_decoupling')
""")
md(r"""
Every plan sets its lamina from its own $D$ ($\delta^* = D/\kappa_g$, so the planned $D$ equals the target by
construction) and meets its $F$ to within one row's extraction; $N D^2$ is constant to about half
a percent at each $F$ (App. B, P-C8.2 reports the same invariance; the alignment matrix records 0.51 % for the cubic
pair over $D$ = 3-7 um, AA23). What a drawing realises is checked by building it: in the release's own builds of
`examples/decoupling/` (`validation/RELEASE.json` lists all nine, each with exit 0 and a passing physics audit) each
drawn mask realises its own $D$ within ±1 % at every gap (Q2, section 6) and its own $F$.
""")

md(r"""
### 4.4 From flows to widths: the closure (App. A Result 10, eq. 69; Result 14)

The march gives flows, a photomask needs widths. Gifford's assumption (pressure equilibration, PE) is that each gap
row equalises the pressure across the device, so in one stretch the side and center segments see the same pressure
drop and carry flows in the ratio of their conductances:

$$ \sigma(i) = \frac{S(i)}{P(i)} = \frac{R_c(i)}{R_s(i)} \qquad\text{(App. A, eq. 69)}. $$

$R$ is the resistance of a rectangular duct (App. A, Result 14: Gifford's eqn 5, the Yang series; the design code's
`r`, opencif's `kernels.r`, and `schedule.rr` for widths in micrometres). At fixed depth it strictly decreases as the width grows (Result 14, "Monotonicity"),
so each $R_s(i)$ gives exactly one side width, found by a one-dimensional search (Gifford increased $w_s$ step by step
until the resistance matched; the cell below uses Brent's method).

opencif does not assume PE: its **exact gap-drop closure** (design specification v5, section 13) keeps the
pressure drop across each gap, and reduces to the PE closure when the gaps' resistance is negligible
(`docs/ALIGNMENT.md` AA17: "GENERALISED"). The cell below plans the straight Gifford-style example
`examples/gifford_physics_linear.yaml` and compares opencif's side widths with the PE closure computed here from
opencif's own $P(i)$ and $S(i)$.
""")
code(r"""
from opencif import spec as SP
from scipy.optimize import brentq
d_lin = SP.load(str(TREE / 'examples' / 'gifford_physics_linear.yaml'))
sp = SP.to_engine(d_lin)[0]
plan = S.plan(sp)
print('plan: N = %d rows, F = %.5f, w_c = %.0f um, h = %.0f um, delta* = %.4f um (D = %.2f um)'
      % (plan.N, plan.F, sp.wc0, sp.h, sp.delta_star, S.KAPPA_G * sp.delta_star))
print('square duct: R h^4 / (mu L) = %.4f (App. A Result 14: 28.4542)' % (S.rr(20.0, 20.0) * (20.0 * S.UM) ** 4 / S.MU))
seg = np.arange(plan.N - 1)                          # segments 1 .. N-1 (between rows); the last is the outlet stretch
w_pe = np.array([brentq(lambda w: S.rr(w, sp.h) - S.rr(plan.wc[i], sp.h) * plan.P[i] / plan.S[i], 1e-3, 5e3) for i in seg])
rel = plan.ws[seg] / w_pe - 1
print('opencif w_s(1), w_s(N-1): %.3f, %.3f um; PE closure: %.3f, %.3f um' % (plan.ws[0], plan.ws[seg[-1]], w_pe[0], w_pe[-1]))
print('largest relative difference, exact gap-drop vs PE closure: %.1e' % np.max(np.abs(rel)))
S.clear_caches()
g_lin = PL.prepare(SP.load(str(TREE / 'examples' / 'gifford_physics_linear.yaml')), lambda m: None)[3]
m13, f_row = g_lin['M13_gap_dominance_max'], plan.f[5]
rs_rel = plan.Rs[seg] / (plan.Rc[seg] * plan.P[seg] / plan.S[seg]) - 1
print('gap dominance M13 = %.4f (gate 0.05); R_s(exact) / R_s(PE) - 1 = %.4e (median over the rows; %.4e to %.4e); '
      '-2 f M13 = %.4e' % (m13, np.median(rs_rel), rs_rel.min(), rs_rel.max(), -2 * f_row * m13))
fig, ax = plt.subplots(1, 2, figsize=(10, 3.2))
ax[0].plot(seg + 1, plan.ws[seg], color='#2a78d6', label='side width $w_s(i)$')
ax[0].axhline(sp.wc0, color='#eb6834', lw=1, label='center width $w_c$')
ax[0].set_xlabel('gap row i'); ax[0].set_ylabel('um'); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
ax[1].plot(seg + 1, 1e6 * rel, color='#2a78d6')
ax[1].set_xlabel('gap row i'); ax[1].set_ylabel('(exact / PE - 1) x 1e6'); ax[1].grid(alpha=0.3)
ax[1].set_title('exact gap-drop closure vs PE closure')
plt.tight_layout(); show(fig, 'fig_4_4_closure')
""")
md(r"""
The two closures agree to a few parts in $10^5$ on this device. That is **not** because the gaps' resistance is
negligible: the gap-dominance number (the design code's M13, $R_{\mathrm{gap}} f / R_c$) is 0.041, close to its 0.05
gate. It is because the exact closure differs from PE only by the *change* of the gap's pressure drop from one row to
the next (spec v5, section 13: $R_c Q_c - R_s Q_s = R_{\mathrm{gap}}Q_{\mathrm{gap}}(i) - R_{\mathrm{gap}}Q_{\mathrm{gap}}(i+1)$).
With a constant $f$ that change is $2f^2$ of the center flow, so $R_s$ moves by $2f\,\mathrm{M13} \approx 4\times10^{-5}$,
as printed above; on a tapered device (exercise 4) $f$ changes from row to row. The first side channel is
only a few micrometres wide: a sealed start asks for a first side width near zero (App. A, Result 10, "Row-1
closure": $R_s(1)\approx R_c/f$). Printing it needs a fine process: this example assumes a 4.5 um floor and a 25:1 aspect-ratio limit (its 5.3 um
side channel is 80 um deep, 15:1; opencif's default limit is 10:1, `docs/ALIGNMENT.md` AG7c), a stated capability,
not a foundry fact. On a 19 um process the generator draws a **meander head** instead (section 6).
""")

md(r"""
## 5. The whole path in 300 lines: a straight ladder by hand

`opencif/teaching.py` does, for the simplest device, everything the generator does: the steps of sections 3 and 4
($\delta^* = D/\kappa_g$, the fraction $f$, the march, the resistances, the closure, the widths), then the drawing,
the photomask, and an audit that reads the photomask back **knowing nothing of the design**. The device is a straight
ladder with a constant center width and a sealed inlet (the example of section 4.4); the module is short enough to
read in one sitting, and each step is named after the result of Appendix A it uses. The cell below runs its design
step and compares it with the generator's plan of the same device.
""")
code(r"""
from opencif import teaching as T
print('opencif/teaching.py: %d lines' % open(TREE / 'opencif' / 'teaching.py').read().count('\n'))
p_ = d_lin['process']
lad = T.design(D_um=d_lin['targets']['cutoff_um'], F=d_lin['targets']['filtrate_fraction'], depth_um=p_['depth_um'],
               gap_um=p_['gap_um'], post_um=p_['post_length_um'], wall_um=p_['post_width_um'],
               centre_um=d_lin['centre_channel']['width_in_um'], outlet_um=sp.L_out)
print('teaching.design: N = %d rows, f = %.4e at every row, w_s from %.3f to %.3f um, F = %.6f'
      % (lad.N, lad.f, lad.ws[0], lad.ws[-1], lad.F))
same = (lad.N == plan.N and np.array_equal(lad.P, plan.P) and np.array_equal(lad.S, plan.S)
        and np.array_equal(lad.Rs, plan.Rs) and np.array_equal(lad.ws, plan.ws[:plan.N - 1]))
print('the same march, closure and side widths as the generator\'s plan, number for number:', same)
""")
md(r"""
Now a smaller ladder ($F$ = 0.3, so about 400 rows) is drawn, written as a photomask, and read back. The audit finds
the posts (the file's ISLANDS layer), the gaps (the openings along a line through the middle of each wall), the center
and side widths (cuts across the channels), computes every resistance (App. A Result 14), solves the drawn device's
flow network (Kirchhoff's laws), and reports each gap's realised cutoff $D_i = \kappa_g\,x(f_i)\,w_c$ and the
filtrate fraction. The depth cannot be seen on a mask, so it is the one input.
""")
code(r"""
small = T.design(D_um=6.0, F=0.3, depth_um=80, gap_um=20.5, post_um=25.3, wall_um=59.5, centre_um=100)
(OUT / 'teaching').mkdir(exist_ok=True)
path_s = OUT / 'teaching' / 'ladder.dxf'
T.write_photomask(T.draw(small), str(path_s))
a_s = T.audit_photomask(str(path_s), depth_um=80, D_um=6.0)
print('the photomask: %d gaps and %d posts per side; read back: wall %.2f um, gaps %.3f-%.3f um, center %.3f-%.3f um'
      % (a_s['gaps_per_side'], a_s['posts_per_side'], a_s['wall_um'], *a_s['gap_um'], *a_s['centre_um']))
print('realised cutoff at every gap: %.6f - %.6f um (D = 6 um; largest relative deviation %.1e)'
      % (*a_s['D_um'], a_s['D_rel_dev_max']))
print('filtrate fraction on the drawn mask %.6f (the design: %.6f)' % (a_s['F'], small.F))
fig, ax = plt.subplots(figsize=(10, 3.4))
fl_ = T.draw(small)
from shapely.geometry import box as _box
win_ = fl_.intersection(_box(-60, -260, 6 * small.pitch, 260))
for g_ in getattr(win_, 'geoms', [win_]):
    xs_, ys_ = g_.exterior.xy
    ax.fill(xs_, ys_, color='#cfe3f7', lw=0)
    for r_ in g_.interiors:
        xi_, yi_ = r_.xy
        ax.fill(xi_, yi_, color='#8a8a8a', lw=0)
ax.set_aspect('equal'); ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)')
ax.set_title('The teaching ladder\'s first rows, drawn from the polygon written to the DXF: fluid (blue), posts (grey)')
show(fig, 'fig_5_teaching_ladder')
""")
md(r"""
**A seeded fault.** The cell below widens one gap of the drawing by 0.6 um (2.9 % of its 20.5 um), writes the
photomask again and audits it blind. The widened gap skims more than its share, so its cutoff rises: the audit finds
it at that row without being told where to look. The physics auditor of section 7 does the same on every photomask
the generator writes, with the full device's parts (its check P5 fails a gap whose lamina exceeds $\delta^*$ by more
than 2 %).
""")
code(r"""
k_ = 200                                                   # gap 201 (0-based 200), side A
x_ = k_ * small.pitch
bad_ = T.draw(small).union(_box(x_ + small.G / 2 - 0.01, small.w_c / 2 + 0.5, x_ + small.G / 2 + 0.6,
                                small.w_c / 2 + small.l_g - 0.5))
path_b = OUT / 'teaching' / 'ladder_one_gap_wider.dxf'
T.write_photomask(bad_, str(path_b))
a_b = T.audit_photomask(str(path_b), depth_um=80, D_um=6.0)
dev_ = a_b['D_rows'][:, 0] / 6.0 - 1
print('largest deviation of the realised cutoff: %+.2f %% at gap %d (side A)' % (100 * np.nanmax(dev_), int(np.nanargmax(dev_)) + 1))
fig, ax = plt.subplots(figsize=(8, 3.0))
ax.axhspan(-1, 1, color='#efeeea', label='Q2 tolerance, +-1 %')
ax.plot(np.arange(1, len(dev_) + 1), 100 * dev_, color='#2a78d6', lw=0.9, label='one gap 0.6 um wider')
ax.plot(np.arange(1, len(dev_) + 1), 100 * (a_s['D_rows'][:, 0] / 6.0 - 1), color='#8a8a8a', lw=0.9, ls='--', label='as designed')
ax.set_xlabel('gap row i'); ax.set_ylabel('$D_i / D - 1$ (%)'); ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.set_title('The blind audit finds the widened gap')
show(fig, 'fig_5_seeded_fault')
""")

md(r"""
## 6. One call: from $(D, F)$ to a photomask

`opencif.generate(D, F)` builds the whole device: it writes the specification it used (every default recorded), plans
the rows, draws the photomask, measures the drawn polygons (every side segment, head unit, turn and outlet lead),
solves the 1-D flow network *on the drawn geometry*, and checks the result. Without other arguments it uses the
`generic19` process (100 um deep, 19 um gaps and floor), one straight leg, a constant center width $w_c = D/0.05$, a
**packed meander head** (long folded channels, 1.45 times the gap wide, that carry the first rows' tiny side flows, so
the start is sealed: gap 1 skims exactly $f^*$, the design fraction), a plain inlet, separate outlets and **strict** mode (`docs/OPTIONS.md`
lists every part and default). Below: $D$ = 4 um, $F$ = 0.5 (about 2-10 minutes, the longer while section 10's builds run beside it).
""")
code(r"""
t0 = time.time()
first = opencif.generate(4.0, 0.5, out=OUT / 'first')
print('%.0f s' % (time.time() - t0))
d = first.record['design']
print('verdict  :', first.verdict)
print('exit code:', first.exit_code, '(0 = PASS)')
print('N = %d gap rows (%d of them in the meander head), F planned %.5f, F on the drawn mask %.5f'
      % (d['N'], d['NM'], d['F'], first.record['validation']['network']['F_as_drawn']))
print('photomask:', os.path.relpath(first.dxf, HERE), '| sha256', sha256(first.dxf)[:16])
print('physics audit of the photomask file:', first.physics_audit)
print('spec written to', os.path.relpath(first.spec_path, HERE), '(rebuild it with: python -m opencif build <that file> --out DIR)')
""")
code(r"""
for fig_path in first.figures:          # the build's own proof figures: the layout, the zoom windows, the physics
    preview(fig_path)
""")
md(r"""
**What the record checks.** The record (`<name>.json`) holds every check; `assessment.status` summarises them. Two
checks carry the physics on the **drawn** mask:

* **Q1**: the network solved on the drawn geometry delivers the designed fraction at every gap;
* **Q2**: the cutoff each gap *realises* on the drawn mask, $D_i$, is within the tolerance (±1 %) of $D$, and the
  delivered $F$ is within its tolerance of the target. In reproduce mode (section 10) a row is exempt from Q2 only when
  a declared departure (a seed gap, a head drawn as depicted, a floor start, the papers' transition, the loop pull-in)
  is shown to move it (`docs/REFERENCE.md`, "What every build checks").

$F$ here is the filtrate fraction of the 1-D network on the drawn mask. 3-D flow at the outlet bifurcation is not
modelled: on Dinh's V1 the measured $F$ was about 10 % below the 1-D value (the glossary's $\eta = 0.903$, not
applied; `docs/LIMITATIONS.md` item 5), and App. B attributes the difference to the outlet's geometry (P-C4.2), so a
fabricated device's $F$ may fall short of the 1-D value by a similar amount (section 11).

The cell below plots $D_i/D - 1$ for both side channels, from the arrays the build saved (`<name>.pkl`; load only
pickles you built yourself).

**The bracket in the verdict is not Q2.** It is the *1-D model band* (`docs/REFERENCE.md`, "What every build checks"): the same
drawn network re-solved with the junction terms (the bends and the head's T-junctions) at the two ends of their
uncertainty, reported and never a pass condition. Here it moves the rows at the seam between the meander head and the
straight section, which are the most sensitive to the head's resistance. Rows "losing their cutoff" in that bracket are rows inside
the tolerance on the drawn mask that, in one of those scenarios, would have reverse flow or skim half the flow or more
($f \ge 1/2$), so that no cutoff is defined there.
""")
code(r"""
import pickle
arr = pickle.load(open(Path(first.out) / (first.record['design']['name'] + '.pkl'), 'rb'))['arrays']
rc = first.record['validation']['realised_cutoff']
dA, dB = np.asarray(arr['dA'], float), np.asarray(arr['dB'], float)
rows_ = np.arange(1, len(dA) + 1)
fig, ax = plt.subplots(figsize=(8, 3.2))
ax.axhspan(-1, 1, color='#efeeea', label='Q2 tolerance, +-1 %')
ax.plot(rows_, 100 * (dA / rc['D_target_um'] - 1), color='#2a78d6', lw=0.8, label='side A')
ax.plot(rows_, 100 * (dB / rc['D_target_um'] - 1), color='#eb6834', lw=0.8, ls='--', label='side B')
ax.axvline(first.record['design']['NM'] + 0.5, color='k', lw=0.6)
ax.text(first.record['design']['NM'] + 8, 0.7, 'end of the meander head', fontsize=8)
ax.set_xlabel('gap row i'); ax.set_ylabel('$D_i / D - 1$ (%)'); ax.set_ylim(-1.5, 1.5); ax.legend(fontsize=8)
ax.set_title('Realised cutoff on the drawn mask'); ax.grid(alpha=0.3)
show(fig, 'fig_5_realised_cutoff')
print('realised cutoff range %.4f - %.4f um (D = %.1f um); share of F through rows within 1 %% of D: %.3f'
      % (rc['D_realised_min_um'], rc['D_realised_max_um'], rc['D_target_um'], rc['share_of_F_within_1%']))
print('status:')
for k, v in first.record['assessment']['status'].items():
    print('  %-30s %s' % (k, v[:110]))
""")
md(r"""
**The photomask itself.** The DXF has two layers: `CHANNELS` (the fluid region, one polygon) and `ISLANDS` (the posts
and walls inside it). Below, the file is read with `ezdxf` and a window of the straight section is drawn from it.

**Exit codes** (`docs/ARCHITECTURE.md`, section 8): 0 PASS (or REPRODUCED); 2 FAIL (kept as `<name>_FAILED.dxf`); 3 outside the validity envelope
in evaluate or reproduce mode; 4 INCOMPLETE; 5 REFUSED (no photomask written; the record `<name>.json` gives the
reasons and the cutoff band); 6 an invalid specification; 1 the generator could not draw the request (an
unanticipated exception; no record).
""")
code(r"""
import ezdxf
from matplotlib.patches import Polygon as MplPolygon
doc = ezdxf.readfile(first.dxf)
msp = doc.modelspace()
layers = {}
for e in msp:
    layers.setdefault(e.dxf.layer, []).append(e)
print({k: len(v) for k, v in layers.items()}, 'entities per layer')
from shapely.geometry import Polygon as ShPolygon
outline = [e for e in layers['CHANNELS'] if e.dxftype() == 'LWPOLYLINE'][0]
fluid = ShPolygon([(p[0], p[1]) for p in outline.get_points('xy')])
x0_, y0_, x1_, y1_ = fluid.bounds
print('the fluid outline: %d vertices; the whole mask %.1f x %.1f mm' % (len(fluid.exterior.coords), (x1_ - x0_) / 1000, (y1_ - y0_) / 1000))
""")
md(r"""
The polygons are in micrometres. Rather than guessing where the rows are on the chip, the cell below uses the drawn
mask's own coordinates: it takes the posts (ISLANDS smaller than a wall) and zooms on the one in the middle of the device.
""")
code(r"""
isl = [ShPolygon([(p[0], p[1]) for p in e.get_points('xy')]) for e in layers.get('ISLANDS', [])
       if e.dxftype() == 'LWPOLYLINE' and len(list(e.get_points('xy'))) > 2]
small = sorted([p for p in isl if p.area < 20000], key=lambda p: p.centroid.x)     # posts, not walls
cx = np.median([p.centroid.x for p in small])
mid = min(small, key=lambda p: abs(p.centroid.x - cx))
win = (mid.centroid.x - 400, mid.centroid.x + 400, mid.centroid.y - 250, mid.centroid.y + 250)
fig, ax = plt.subplots(figsize=(9, 5.6))
for lay, col in (('CHANNELS', '#cfe3f7'), ('ISLANDS', '#8a8a8a')):
    for e in layers.get(lay, []):
        if e.dxftype() != 'LWPOLYLINE':
            continue
        pts = np.array([(p[0], p[1]) for p in e.get_points('xy')])
        if len(pts) > 2 and pts[:, 0].max() > win[0] and pts[:, 0].min() < win[1] and pts[:, 1].max() > win[2] and pts[:, 1].min() < win[3]:
            ax.add_patch(MplPolygon(pts, closed=True, color=col, lw=0))
ax.set_xlim(win[0], win[1]); ax.set_ylim(win[2], win[3]); ax.set_aspect('equal')
ax.set_xlabel('x (um)'); ax.set_ylabel('y (um)')
ax.set_title('The photomask, drawn from the DXF: posts (grey) and fluid (blue), 800 x 500 um window')
show(fig, 'fig_5_photomask_window')
""")

md(r"""
## 7. The physics audit: the photomask read back from its file

The checks of section 6 run on the polygon the generator drew. The **physics audit** then reads the photomask back
from its DXF file and checks the physics on what the file contains. It is the design code's own physics auditor
(`CIF_v2.py`'s embedded `ladder_physics.audit`, checks P1-P10) ported to opencif (`opencif/checks/auditor.py`;
`docs/PHYSICS_AUDIT.md` gives each check's source and what the port changes): Appendix B says that the C7/C8 software
is "warranted by their acceptance suite and the physics auditor" and that "the physics auditor separately checks the
drawn device's inlet against its own, stricter channel limit" (App. B, Scope).

| check | what it requires |
|---|---|
| P0 | the file read back is the drawn polygon, vertex for vertex |
| P1 | the kernels reproduce their published anchors |
| P2 | mass is conserved: in the plan, and at every node of the drawn network |
| P3 | the drawn side channels realise the closure |
| P4 | in the bulk of the device every gap skims the kernel's fraction, within 0.5 % |
| P5 | no gap's lamina exceeds $\delta^*$ by more than 2 %: the cutoff guarantee |
| P6, P6b | the filtrate fraction: on target, and equal to App. A's closed form at the outlet |
| P7 | how much the outlet's matching matters (reported) |
| P8 | the auditor's own, stricter inlet limit, Re $\le$ 1 (reported: the design envelope is Re $\le$ 10) |
| P8-M14 | the envelope's Reynolds and Stokes gates, on the drawn mask |
| P9, P9-AR | the fabrication floors and the aspect ratio |
| P10 | no corner sharper than the process allows |

The first device's audit, from its record:
""")
code(r"""
pa = first.record['validation']['physics_audit']
print('physics audit of %s: %s (%.1f s)' % (pa['photomask'].get('file_final', pa['photomask']['file']), pa['status'], pa['seconds']))
for c in pa['checks']:
    v = c['value']
    short = {k: v[k] for k in list(v)[:3]} if isinstance(v, dict) else v
    print('%-7s %-9s %-8s %s' % (c['id'], 'ok' if c['ok'] else ('reported' if c['role'] == 'report' else 'FAIL'),
                                 c['role'], json.dumps(short, default=float)[:110]))
print('P8, the auditor\'s own limit:', next(c for c in pa['checks'] if c['id'] == 'P8')['value']['result'])
""")
md(r"""
P8 says what App. B's sentence says: the auditor checks the drawn inlet against its own, stricter limit (Re $\le$ 1).
At this device's flow, 30 uL/min, its inlet runs at Re 5.6: inside the design envelope's channel gate (Re $\le$ 10,
P8-M14) and above the auditor's own limit, which it would meet at 5 uL/min (the auditor's default flow). The audit
records both; only the envelope's gate decides the verdict (`docs/PHYSICS_AUDIT.md`, `docs/LIMITATIONS.md` item 6).

**The audit from scratch.** `python -m opencif audit <folder>` (or `auditor.audit_run`) repeats the audit of a
finished build starting from the record's input alone: it plans and draws the design again, reads the DXF, requires
it to equal that drawing, measures the polygon **read from the file**, runs the DRC on the file, and loads no saved
Python object. It must give the same value for every check (about a minute here).
""")
code(r"""
from opencif.checks import auditor as AU
t0 = time.time()
again = AU.audit_run(first.out, say=lambda m: None)
same_ = all(a_['id'] == b_['id'] and a_['ok'] == b_['ok'] and a_['value'] == b_['value']
            for a_, b_ in zip(pa['checks'], again['checks'])) and len(pa['checks']) == len(again['checks'])
print('standalone audit: %s in %.0f s; every check\'s value equal to the build\'s own audit: %s'
      % (again['status'], time.time() - t0, same_))
""")
md(r"""
`tests/test_auditor.py` runs the port beside the design code's own auditor on the design code's demonstration
devices (the same numbers), and seeds faults into a drawn device one at a time: a vertex of the file moved by 1 nm,
a kernel off by $10^{-9}$, a gap widened by 0.6 um, a side segment narrowed, the flow tripled, a lost feature, a sharp
corner, a raised target $F$, a tighter aspect-ratio limit. Each is caught by the check meant to catch it.
""")

md(r"""
## 8. The validity envelope: what the generator refuses, and why

The cutoff law and the 1-D model hold only inside an **envelope** (App. B, thesis and scope; the design code's gates;
`docs/REFERENCE.md`, "What every build checks"). In **strict** mode (the default) a design outside it is **refused before
anything is drawn** (exit 5), with the reason. The gates include:

* the $\kappa$ plateau, $\tilde d = D/w_c \ge 0.035$ at the widest center width, and the foot (App. B C3c);
* the channel Reynolds number $\mathrm{Re} \le 10$ and the Stokes gates $\mathrm{Re}_{\mathrm{gap}}\le10^{-2}$,
  $\mathrm{Re}_\delta\le10^{-2}$ (the Reynolds numbers of the flow through a gap and of the skimmed lamina; App. B's thesis: "the Stokes regime is
  gated in code"; the design code's M14);
* gap dominance $R_{\mathrm{gap}}f/R_c < 0.05$ (the design code's M13, spec v5 section 19), the thin lamina
  $\delta/h \le 0.10$ (E11), a quasi-developed taper (M16), and the process floor
  $w_{\min}=\max(F_{\min}, h/AR_{\max})$ on the gap, the head's channel, every center width and every straight side
  channel (G7); on the drawn mask the narrowest channel's aspect ratio is then a pass condition (exit 2 in strict
  mode).

In **evaluate** mode the same design is drawn and reported (exit 3); **reproduce** mode is evaluate mode plus the
papers' deliberate departures (a seed gap, a head drawn as depicted, a floor start, the papers' transition, the loop
pull-in), each
declared row by row (section 10). The cells below *plan* three refused designs (seconds each).
""")
code(r"""
cases = [('off the plateau: D 2 um on a 100 um center channel', dict(cutoff=2.0, F=0.5, taper=('constant', 100))),
         ('too fast: 200 uL/min', dict(cutoff=4.0, F=0.5, flow_uL_min=200)),
         ('gap dominance and Stokes: D 6 um on 80 um, 100 um deep', dict(cutoff=6.0, F=0.5, taper=('constant', 80)))]
for k, (title, kw) in enumerate(cases, 1):
    kw = dict(kw); Dv, Fv = kw.pop('cutoff'), kw.pop('F')
    parts = dict(DC); parts.update(kw)
    r = opencif.generate(Dv, Fv, out=OUT / 'plans', check_only='plan', name='refused_%d' % k, **parts)
    print('%s\n   exit %d: %s\n' % (title, r.exit_code, r.verdict))
""")
code(r"""
import yaml
# the same straight Gifford-style example at 30 uL/min: strict mode refuses it (the lamina's Reynolds number at gap 1 is
# just above 1e-2); evaluate mode plans it and reports it (exit 3)
for mode in ('strict', 'evaluate'):
    d_ = yaml.safe_load(open(TREE / 'examples' / 'gifford_physics_linear.yaml'))
    d_['operation']['flow_uL_min'] = 30
    d_['gates']['mode'] = mode
    d_['name'] = 'gifford_linear_30_' + mode
    r = opencif.generate(d_['targets']['cutoff_um'], d_['targets']['filtrate_fraction'], process=d_['process'],
                         chip=d_['chip'], taper=d_['centre_channel'], layout=d_['layout'], meander=None,
                         inlet=d_['inlet'], outlet=d_['outlet'], flow_uL_min=30, mode=mode, name=d_['name'],
                         out=OUT / 'plans', check_only='plan')
    print('%-8s exit %d: %s' % (mode, r.exit_code, r.verdict))
""")

md(r"""
## 9. Composing new devices from parts

A device is composed of parts, each a section of the specification (`docs/OPTIONS.md` gives each part's styles,
origin, an example and the record field that verifies it on the drawn mask):

| part | styles | origin |
|---|---|---|
| `layout` | `linear`, `serpentine` (legs side by side), `paperclip` (nested legs) | Gifford 2014 Fig. 1 (linear); Gifford 2014 Fig. 4 and Xia 2016 Fig. 3a (paperclip); this project (serpentine) |
| `centre_channel.taper` | `constant`, `constant_shear` (narrows as the flow drops, holding the wall shear) | Gifford 2014; Xia 2016 |
| `meander` (the head) | `pack`, `arms`, `hairpins`, `narrow`, or none | this project; Dinh 2024; Xia 2016; this project |
| `inlet` | `plain`, `filter`, `dinh_filter` | this project (filter); Dinh 2024 Fig. 2A (dinh_filter) |
| `outlet` | `separate`, `teardrop`, `ring` | Gifford 2014 Fig. 4; Dinh 2024 Fig. 2A; after Xia 2016 Fig. 3a |
| `start` | `sealed`, `floor` (no head: the first side channel at the process floor) | Gifford 2014 Fig. 2 (sealed) and p. 4498 (floor); App. A Result 10 |
| `array` | n copies of a passing element, one inlet | this project (`opencif.compose`) |

Every part is resistance-exact and checked gap by gap on the drawn mask. Three new devices follow: a serpentine
platelet-scale device on a microscope-slide-sized chip, a nested paperclip with a ring outlet (one of the release's
own test designs, rebuilt here and compared with the release's photomask byte for byte), and an array of three copies.
""")
code(r"""
t0 = time.time()
serp = opencif.generate(2.5, 0.6, process='dinh2024', chip=(25, 75), layout='serpentine',
                        taper=('constant_shear', 70, 50), meander='pack', inlet='filter', outlet='teardrop',
                        out=OUT / 'serpentine', name='serp_D2p5_F0p6')
print('%.0f s | exit %d | %s' % (time.time() - t0, serp.exit_code, serp.verdict))
dd = serp.record['design']
print('N = %d rows (%d in the head), %d legs, F on the drawn mask %.5f; over the chip area: %s'
      % (dd['N'], dd['NM'], dd.get('n_legs', 1), serp.record['validation']['network']['F_as_drawn'], serp.over_chip_area))
preview(serp.figures[0])
""")
code(r"""
name = 'm11_paperclip_const_pack_filter_ring'
spec_path = TREE / 'examples' / 'matrix' / (name + '.yaml')
print(open(spec_path).read())
t0 = time.time()
rec = PL.run(str(spec_path), str(OUT / name), verbose=False)
dxf = OUT / name / (name + '.dxf')
print('%.0f s | exit %d | %s' % (time.time() - t0, rec['assessment']['exit_code'], rec['assessment']['verdict']))
print('photomask sha256 %s\nreleased         %s\nbyte-identical to the release\'s own build: %s'
      % (sha256(dxf), FROZEN[name]['dxf_sha256'], sha256(dxf) == FROZEN[name]['dxf_sha256']))
print('physics audit:', rec['validation']['physics_audit']['status'])
preview(OUT / name / (name + '_layout.png'))
""")
code(r"""
# an array: three copies of a passing element (D 5 um, F 0.5, straight, packed head) sharing one inlet; each copy keeps
# its own outlets, so each realises the element's D and F whatever the inlet manifold does (docs/REFERENCE.md, Composition)
t0 = time.time()
arr3 = opencif.generate(5.0, 0.5, array=dict(count=3, inlet='shared'), out=OUT / 'array', name='arr_D5_F0p5')
print('%.0f s | element: exit %d | %s' % (time.time() - t0, arr3.exit_code, arr3.verdict))
a_ = arr3.record.get('array') or {}
print('array:', a_.get('verdict', a_))
for f_ in sorted(glob.glob(str(OUT / 'array' / '*_x3*.png'))):
    preview(f_)
""")

md(r"""
## 10. Reproducing the published devices

The release carries **presets** for devices of the three papers. Every value in a preset carries its source
(`docs/PRESET_PROVENANCE.md`): **P** printed, **M** measured on a figure, **D** derived through the frozen model or
fitted to a measurement, **A** assumed, **G** the generator's choice. The reproductions are **design-level**: they
reproduce the papers' design quantities and align with their micrographs; no device was fabricated or measured
(`docs/REPRODUCTIONS.md`). Every preset runs in **reproduce** mode: the device is drawn with the paper's deliberate
departures, each declared and attributed row by row (`docs/OPTIONS.md`, `gates.mode`). A preset outside the validity
envelope is drawn and reported with exit 3: here Gifford's Fig. 6 ($\tilde d_{\min}$ 0.029) and Xia's RS-7
($\tilde d_{\min}$ 0.021; also M13, the channel Reynolds number and the Stokes gates). Dinh's figure-aligned preset (its cutoff fitted to
3.67 um) and Gifford's Fig. 4 lie inside it and reproduce with exit 0. At 2.17 um, the cutoff that Dinh's printed $f^*$ = 1.04 × 10⁻⁴ at $w_c$ = 75 um implies through
$\kappa_g$, the plateau guard would refuse the device, as it would Xia's RS-7 (App. B, C8 CQ7: Dinh's widest end "gives u = d/w_c = 0.029 ...
and the plateau guard refuses it").

| preset | paper | what is reproduced |
|---|---|---|
| `dinh2024_v1_figure_aligned` | Dinh 2024, V1 (platelets) | the printed design, $D$ fitted to the printed length |
| `gifford2014_fig4_aligned` | Gifford 2014, Fig. 4 | the nested four-leg device of the inset (register GIF4-1: leg 2's widths remain open) |
| `gifford2014_fig6_aligned` | Gifford 2014, Fig. 6 (platelets) | seven nested legs; the start fitted on leg 1 |
| `xia2016_rs7_aligned` | Xia 2016, RS-7 (leukocytes) | the constant-shear taper and hairpin head, $D$ and $F$ fitted to Fig. 3 |
""")
code(r"""
head = open(TREE / 'presets' / 'gifford2014_fig4_aligned.yaml').read().split('\nname:')[0]
print(head[:2600], '\n...')
""")
code(r"""
import yaml
for n in REPRO:
    t0 = time.time()
    dsp = yaml.safe_load(open(TREE / 'presets' / (n + '.yaml')))
    S.clear_caches()
    spec_, opts_, sc_, g_, env_bad_, notes_ = PL.prepare(SP.load(dsp), lambda m: None)
    print('%-28s planned in %4.1f s: N %5d, head rows %3d, F %.5f, D %.4f um; outside the envelope: %s'
          % (n, time.time() - t0, sc_.N, sc_.NM, sc_.F, S.KAPPA_G * spec_.delta_star, '; '.join(env_bad_) or 'no'))
    print('%-28s the release\'s build: exit %d, N %d, F %.5f, physics audit %s' % ('', FROZEN[n]['exit_code'], FROZEN[n]['N'],
                                                                                FROZEN[n]['F_design'], FROZEN[n]['physics_audit']))
""")
md(r"""
**The full builds**, started in section 1: each preset is built in its own process
(`python -m opencif build presets/<name>.yaml`), two at a time; a build already finished by an earlier run of this
notebook is not repeated if its photomask's sha256 is the released one (section 14 then records which builds were
reused, and the code stamp they were made with). The cell below waits for them; the next one
compares each photomask with the one the release's own rebuild wrote: they must be identical byte for byte, and
each must pass its physics audit (a reproduction outside the validity envelope passes it with P8-M14 outside, exit 3).
""")
code(r"""
if repro_thread is not None:
    t_wait = time.time()
    repro_thread.join()
    print('waited %.0f min for the background builds' % ((time.time() - t_wait) / 60))
print(REPRO_LOG.read_text() if REPRO_LOG.exists() else '(no builds: quick mode)')
""")
code(r"""
print('%-28s %-5s %-6s %-9s %-18s %-22s %s' % ('preset', 'exit', 'N', 'F drawn', 'photomask sha256', 'physics audit', 'equal to the release\'s'))
for n in REPRO:
    rj = OUT / 'repro' / n / (n + '.json')
    if not rj.exists():
        print('%-28s (not built: quick mode)' % n)
        continue
    r_ = json.load(open(rj))
    dxf = OUT / 'repro' / n / r_['files']['dxf']
    h = sha256(dxf)
    print('%-28s %-5d %-6d %-9.5f %-18s %-22s %s' % (n, r_['assessment']['exit_code'], r_['design']['N'],
          r_['validation']['network']['F_as_drawn'], h[:16], r_['validation']['physics_audit']['status'],
          'yes, byte for byte' if h == FROZEN[n]['dxf_sha256'] else 'NO'))
    print('    verdict: %s' % r_['assessment']['verdict'][:240])
""")
md(r"""
**The dissertation figures** for each reproduction (`tools/dissertation_figures.py`, outside the code stamp: it draws
what a build holds and computes nothing a record does not): the photomask with four zoom windows; the physics panel
(the realised cutoff per gap with the ±1 % band and the declared rows; the cumulative filtrate fraction; the 1-D
model's scenarios); a provenance table. None contains a paper's image.

The provenance tables give the realised cutoff over every gap, declared rows included. Dinh's figure-aligned preset
reproduces (exit 0) although its realised cutoffs span 2.2-28.7 um: all 232 of its rows outside ±1 % are declared
(gap 1's seed, the 217 rows of the head drawn as depicted, and the 14 rows of the papers' transition, six of them with
reverse flow), and none is undeclared. Xia's realised minimum, 0 um, is a row of its head drawn as depicted; its 2476
rows outside ±1 % are likewise all declared (most of them the loop pull-in).
""")
code(r"""
FIGS = OUT / 'dissertation_figures'
FIGS.mkdir(exist_ok=True)
for n in REPRO:
    run_dir = OUT / 'repro' / n
    if not (run_dir / (n + '.json')).exists():
        continue
    p = subprocess.run([sys.executable, '-I', str(TREE / 'tools' / 'dissertation_figures.py'), 'build', str(run_dir), str(FIGS)],
                       capture_output=True, text=True)
    said = (p.stdout.strip().splitlines() or [''])[-1] if p.stdout else p.stderr[-300:]
    print(n, 'figures: exit', p.returncode, said.replace(str(HERE) + os.sep, ''))
for n in REPRO:
    for kind in ('overview', 'physics'):
        f_ = FIGS / ('%s_%s.png' % (n, kind))
        if f_.exists():
            display(Markdown('**%s, %s**' % (n, kind)))
            preview(f_, width=1000)
    pv = FIGS / ('%s_provenance.md' % n)
    if pv.exists():
        display(Markdown(pv.read_text()))
""")
md(r"""
**The taper-direction test (V13) on Gifford's Fig. 4.** In a CIF leg the side channels widen along the flow, so the
sign of a leg's outer width slope is its flow direction (Gifford 2014, Fig. 1 caption). RC8 measured that slope on the
Fig. 4 inset for the neighbouring leg, through each candidate layout's registration, and compares it with the
photomask's: the nested layout (legs 2 and 4, as the caption says) passes, the path-order layout fails. The figure
below is drawn from RC8's frozen trace records (widths only: no paper image).
""")
code(r"""
G4 = TREE / 'provenance' / 'rc7_traces' / 'results' / 'gifford'
p = subprocess.run([sys.executable, '-I', str(TREE / 'tools' / 'dissertation_figures.py'), 'v13', str(FIGS / 'fig_taper_direction'),
                    '%s=reading B (path order): the neighbour, leg 3' % (G4 / 'v13_pathorder_neighbour.json'),
                    '%s=the nested reading (the device): the neighbour, leg 2' % (G4 / 'v13_aligned_neighbour.json'),
                    '%s=the control: leg 4\'s array part' % (G4 / 'v13_aligned_leg4.json')], capture_output=True, text=True)
print('exit', p.returncode)
preview(FIGS / 'fig_taper_direction.png', width=1100)
for k in ('aligned_neighbour', 'pathorder_neighbour', 'aligned_leg4'):
    v = json.load(open(G4 / ('v13_%s.json' % k)))
    print('%-20s %s: image slope %+.5f +- %.5f um/um (%.1f SE), mask %+.5f, %d of %d eligible positions, window +-%g um'
          % (k, v['verdict'], v['slope_image'], v['se_image'], abs(v['z_image']), v['slope_mask'], v['positions'],
             v['eligible'], v['window_um']))
""")
md(r"""
**What the reproductions do not show** (the register `RC8_DISCREPANCIES.md`, in the release archive; `docs/LIMITATIONS.md` item 28 for GIF4-1): on Gifford's
Fig. 4 the neighbouring leg's side channels read about 18 um per side wider on the image than the printed schedule
draws them (GIF4-1, open); on Xia's Fig. 3 panel f the side-channel walls lie 36-40 um inside the aligned mask's
(INV-4, open). Both are stated, not hidden; neither changes the generator.
""")

md(r"""
## 11. opencif and the design code `CIF_v2.py`

Appendix A's table of results names the program behind Results 10-14: "the design code `CIF_v2.py`". It is inside
the release, byte-frozen (`opencif/reference/CIF_v2.py`, md5 `0e0961a5...`), and opencif imports its kernels: the
duct resistance `r` and its inverse, the cubic kernel `g_cubic` and its inverse `x_of_f`, the $\kappa$ pairing, the
shear coefficient, the entrance length, the stadium post shape, and its process and tolerance records. The cell below
runs the design code's own reference device: $D$ = 2.17 um (the cutoff that Dinh's printed $f^*$ = 1.04 × 10⁻⁴ at $w_c$ =
75 um implies through $\kappa_g$, at leading order: the glossary's 2.170; Dinh prints no cutoff in micrometres) and $F$ = 0.5592 (App. B C4's
1-D prediction $F_{\mathrm{closed}}$ for Dinh's V1; Dinh measured 0.505). Its center channel narrows from 60 to 42 um,
Dinh's 75 to 52 um rescaled, so that it passes the plateau guard that refuses Dinh's own widths
(`CIF_v2.lineage_request`; section 10).

How the two programs relate (`docs/ALIGNMENT.md`; `foundation/README.md`, open item D95):

| | `CIF_v2.py` (the design code) | opencif RC8 (the generator) |
|---|---|---|
| kernels ($g$, $x(f)$, $r$, $\kappa_g$) | defines them | imports them from `CIF_v2.py` |
| the march (Result 10) | `march` | `cifsched.build` (the trace suite checks App. A's example) |
| widths from flows | `design()`: the PE closure (spec v5 section 9; App. A eq. 69). Its embedded physics core (`compile_unified`): an exact flow-prescribed closure, a Kirchhoff solve of the compiled pattern and a width-correction loop | the exact gap-drop closure, which reduces to eq. 69 (AA17; section 4.4); opencif does not use the design code's physics core |
| the outlet | PE outlet identity (Result 8) | exact outlet closure (AA8, AA20) |
| geometry | linear, serpentine and paperclip layouts with U-turns, folded side channels, stadium posts, ports and leads (`Layout`, `build_mask`); a meander-head compiler (`compile_meander`, `write_dxf_meander`: heads of arms, taper and decay regions, filleted, with 19 geometry validators) | its own heads (the papers' arms and hairpins, the packed and narrow heads), filter inlets, teardrop and ring outlets, arrays; every bend and junction corrected by exact solves; all measured on the drawn mask |
| checks | validator battery K/E/C/P/M/D/G/N/S/R and further families for its embedded cores and staged pipeline; mask checks G17-G23 on its own builds; physics auditor P1-P10 | acceptance suite (12 test suites); every drawn mask measured (Q1, Q2, DRC); P8's criteria reported, the auditor not run on the mask (LIMITATIONS 6) |
| realised $F$ and $D$ | reports $F_{\mathrm{real}} \approx \eta F$ ($\eta$ = 0.903, measured on Dinh's V1) and a depth transfer bias on $D$, never applied | neither: $F$ is the 1-D network's fraction on the drawn mask (LIMITATIONS 5; ALIGNMENT AG7i); opencif neither reports nor applies the transfer bias |
| published devices | Gifford's Fig. 2 design curves (`gifford_fig2`); the lineage geometry, refused by its plateau guard (`lineage_request`) | presets of eight published designs; four built and traced on the papers' micrographs |
""")
code(r"""
t0 = time.time()
dev = C.design(C.reference_request())
print('CIF_v2 %s, design() in %.1f s' % (C.__version__, time.time() - t0))
print(dev.report())
""")
md(r"""
**Two lines of this report have no counterpart in opencif's records.** $\eta = 0.903$ is the measured ratio of
Dinh's filtrate fraction (0.505) to the 1-D prediction (0.5592; App. B, P-C4.2: 3-D flow at the outlet bifurcation
closes the difference): opencif's $F$ is the 1-D value; the difference depends on the outlet's geometry, so a fabricated device's $F$ may
fall short of it by a similar amount (`docs/LIMITATIONS.md` item 5). The transfer bias is the design code's estimate of how the cutoff changes when $\kappa$, fitted on the calibration
array's channel shape, is used at this device's aspect ratio: by that estimate the true cutoff would be about 15 %
above the predicted one. The design code reports it and never applies it; opencif neither reports nor applies it.
Appendix B keeps the depth effect as a bounded factor folded into $\kappa_{\mathrm{app}}$ (0.81-0.98), and with the
cubic kernel the calibration depth does not enter the row count (glossary; `docs/ALIGNMENT.md`, section 5).
""")

md(r"""
## 12. Micrograph traces (the companion notebook)

opencif RC8 registers its photomasks on the papers' own micrographs (Gifford Fig. 4 and Fig. 6, Xia Fig. 3, Dinh
Fig. 2B) and measures how well they align. Those figures are the publishers': Gifford 2014 and Dinh 2024 are
© The Royal Society of Chemistry (no Creative Commons licence); Xia 2016 is CC BY 4.0. They are not in
`opencif_rc8.py` and not in this notebook. The companion notebook `opencif_rc8_traces.ipynb` takes the figures from a
folder you supply (checked by sha256), assembles each device package from the release, and runs its build, trace,
measurement and rendering, comparing every value with RC8's frozen records. Its outputs show the papers' images:
dissertation evidence only, not for a public repository.
""")

md(r"""
## 13. Exercises

1. **Decoupling.** At $F$ = 0.5 and the parts of section 4.3, predict $N$ for $D$ = 3.5 um from $N D^2$; check it with
   a plan. Then change $F$ to 0.6 at $D$ = 3.5 um: which of $N$ and the planned cutoff moves?
2. **The plateau.** For $D$ = 3 um, what is the widest center channel on the plateau? Plan a design 2 um wider in
   strict and in evaluate mode and read the verdicts.
3. **Why a head?** Plan $D$ = 4 um, $F$ = 0.5 with `meander=None` on the `generic19` process. Read the error it raises
   (a `SpecError`, exit 6 on the command line). Which width falls below which floor, and what remedies does the message
   list?
4. **The closure.** Repeat the comparison of section 4.4 for `examples/matrix/m15_linear_shear_nomeander_filter_teardrop.yaml`
   (a tapered center channel). Is the PE closure still a good approximation? Look at the design code's M13 number
   in the plan's gates (`PL.prepare`).
5. **Compose.** Build a paperclip with a ring outlet for $D$ = 5 um, $F$ = 0.7 on a 25 x 75 mm chip. How many legs
   does the generator choose? Inspect the realised cutoff near the U-turns.
6. **Reproduce.** Open `presets/xia2016_rs7_aligned.yaml`. Which values are printed (P), measured (M), derived or
   fitted (D), assumed (A), or the generator's choice (G)? Why is its exit code 3?
""")

md(r"""
## 14. The reproducibility record

Everything this notebook wrote, with its sha256, and the environment it ran in.
""")
code(r"""
lines = ['opencif %s, code stamp %s, tests %s' % (opencif.__version__, PL.CODE['sha256'], env['tests_sha256']),
         'Python %s on %s; installed = tested versions: %s' % (env['python'], env['platform'], env['matches_lock']),
         'repository: %s, released code stamp %s (%s)' % (TREE, RELEASE['code_stamp'],
                                                           'equal' if RELEASE['code_stamp'] == PL.CODE['sha256'] else 'DIFFERENT')]
if QUICK:
    lines.append('quick mode: the reproductions were planned, not built')
elif REUSED:
    lines += ['reproduction built by an earlier run of this notebook (photomask sha256 verified): %s, code stamp %s, written %s'
              % (n, v['code_stamp'], v['written']) for n, v in REUSED.items()]
else:
    lines.append('every file below was written by this run')
lines.append('')
files = sorted(p for p in OUT.rglob('*') if p.is_file() and p.name != 'MANIFEST.txt')
lines += ['%s  %s' % (sha256(p), p.relative_to(OUT)) for p in files]
(OUT / 'MANIFEST.txt').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines[:lines.index('')]))
print('%d files listed in %s' % (len(files), os.path.relpath(OUT / 'MANIFEST.txt', HERE)))
""")

nb = nbf.v4.new_notebook()
nb['cells'] = cells
nb['metadata'] = {'kernelspec': {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'},
                  'language_info': {'name': 'python'}}
nbf.write(nb, sys.argv[1])
print('wrote %s: %d cells' % (sys.argv[1], len(cells)))
