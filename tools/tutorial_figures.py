"""tutorial_figures.py: draw the explanatory sketch of skimming at one gap for the README and the tutorial.

Usage: python tools/tutorial_figures.py [OUT_DIR]   (default docs/figures)

The sketch has two panels. Panel (a) is a plan view of one gap in the upper wall of the center channel, with the
lamina that leaves through the gap, the dividing streamline that bounds it, and two particles that start against the
wall. Panel (b) is the parabolic velocity profile across the center channel with the lamina's share of the flow
shaded, which is the content of Appendix A, Result 1 (f = 3x^2 - 2x^3 with x = delta/w_c). The geometry is not to
scale, because a real lamina is about one hundredth of the channel's width, but every position that the physics fixes
is computed: the lamina's width follows from x, and the streamlines downstream of the gap are placed by conserving the
flow between each streamline and the wall (the march of Result 10 applied to one gap). The drawing uses matplotlib only
and contains nothing taken from a paper.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, Polygon

# The sketch's proportions, in units of the center channel's half-width (the center line is y = 0 and the upper wall
# is y = 1). The lamina is exaggerated about eight times so that it can be seen.
X_LAM = 0.08                         # x = delta / w_c, the lamina's width as a fraction of the channel's width
DELTA = 2.0 * X_LAM                  # the lamina's width in half-widths, because w_c is two half-widths
F_GAP = 3 * X_LAM ** 2 - 2 * X_LAM ** 3   # the fraction of the center flow that the gap takes (Result 1)
WALL = 0.22                          # the wall's thickness, which is the length of the gap's passage
W_SIDE = 0.46                        # the side channel's width
G0, G1 = 0.0, 0.9                    # the gap's upstream and downstream edges along the flow
XL, XR = -3.1, 2.75                  # the drawn stretch of channel
A_SMALL, A_LARGE = 0.10, 0.30        # particle radii, one smaller and one larger than the lamina
X_SMALL, X_LARGE = -2.35, -1.15      # where the two particles are drawn, upstream of the gap

# Colors: the lamina and everything that leaves with it in orange, what stays in the center channel in blue.
C_FLUID, C_WALL, C_LAM = '#eef3f8', '#8a8f96', '#e69f00'
C_OUT, C_STAY, C_INK, C_MUTED = '#c46f00', '#0b5394', '#222222', '#5f6b78'


def cum_flow(s):
    """Return the share of the center flow between the upper wall and a streamline at distance s (in w_c)."""
    return 3 * s ** 2 - 2 * s ** 3


def s_of_cum(phi):
    """Invert cum_flow on 0 <= s <= 1/2 by bisection, which is ample for a drawing."""
    lo, hi = 0.0, 0.5
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if cum_flow(mid) < phi else (lo, mid)
    return 0.5 * (lo + hi)


def y_after_gap(y0):
    """Return where the streamline that runs at height y0 upstream of the gap runs downstream of it.

    The flow between the streamline and the upper wall loses the gap's fraction, and the center flow loses the
    fractions of both walls' gaps, so the streamline moves toward the wall by an amount that is largest near the wall.
    """
    s0 = (1.0 - y0) / 2.0                                    # distance from the wall in units of w_c
    phi1 = (cum_flow(s0) - F_GAP) / (1.0 - 2.0 * F_GAP)
    return 1.0 - 2.0 * s_of_cum(phi1)


def smooth(t):
    """Return a smooth step from 0 to 1 on 0 <= t <= 1 (zero slope at both ends)."""
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3 - 2 * t)


def bezier(p0, p1, p2, p3, n=120):
    t = np.linspace(0.0, 1.0, n)[:, None]
    p = [np.asarray(q, float) for q in (p0, p1, p2, p3)]
    return (1 - t) ** 3 * p[0] + 3 * (1 - t) ** 2 * t * p[1] + 3 * (1 - t) * t ** 2 * p[2] + t ** 3 * p[3]


def lamina_streamline(k):
    """Return a streamline of the lamina; k = 0 runs along the wall and k = 1 is the dividing streamline.

    Fluid next to the wall turns into the gap near its upstream edge, and fluid near the dividing streamline turns
    near its downstream edge, so the streamlines never cross. In the side channel the new filtrate runs along the
    wall, beneath the filtrate that arrived from the gaps upstream.
    """
    y0 = 1.0 - k * DELTA
    x_turn = G0 + (0.12 + 0.70 * k) * (G1 - G0)              # where the streamline crosses the wall's inner face
    y_side = 1.0 + WALL + 0.05 + (1.0 - k) * 0.13           # its height in the side channel
    head = np.column_stack([np.linspace(XL, x_turn - 0.9, 50), np.full(50, y0)])
    rise = bezier((x_turn - 0.9, y0), (x_turn - 0.3, y0), (x_turn - 0.05, 0.92 + 0.05 * (1 - k)), (x_turn, 1.0))
    turn = bezier((x_turn, 1.0), (x_turn + 0.03, 1.0 + 0.6 * WALL), (x_turn + 0.05, y_side),
                  (x_turn + 0.55, y_side))
    tail = np.column_stack([np.linspace(x_turn + 0.55, XR - 0.15, 30), np.full(30, y_side)])
    return np.vstack([head, rise, turn, tail])


def center_streamline(y0):
    """Return a center-channel streamline that moves toward the wall across the gap by flow conservation."""
    x = np.linspace(XL, XR - 0.15, 400)
    y1 = y_after_gap(y0)
    return np.column_stack([x, y0 + (y1 - y0) * smooth((x - (G0 - 0.4)) / (G1 + 0.8 - (G0 - 0.4)))])


def dividing_streamline():
    """Return the dividing streamline, which starts at the lamina's edge and ends on the downstream post's corner."""
    x = np.linspace(XL, G1, 300)
    return np.column_stack([x, 1.0 - DELTA + DELTA * smooth((x - (G0 - 0.9)) / (G1 - (G0 - 0.9))) ** 1.6])


def arrow_on(ax, path, frac, color, lw=1.2):
    i = int(frac * (len(path) - 2))
    ax.add_patch(FancyArrowPatch(path[i], path[i + 1], arrowstyle='-|>', mutation_scale=11, color=color, lw=lw))


def draw_plan(ax):
    ax.add_patch(Rectangle((XL, 0.0), XR - XL, 1.0, color=C_FLUID, zorder=0))
    ax.add_patch(Rectangle((XL, 1.0 + WALL), XR - XL, W_SIDE, color=C_FLUID, zorder=0))
    ax.add_patch(Rectangle((G0, 1.0), G1 - G0, WALL, color=C_FLUID, zorder=0))
    for x0, x1 in ((XL, G0), (G1, XR)):                              # the posts on either side of the gap
        ax.add_patch(Rectangle((x0, 1.0), x1 - x0, WALL, color=C_WALL, zorder=3))
    ax.add_patch(Rectangle((XL, 1.0 + WALL + W_SIDE), XR - XL, 0.10, color=C_WALL, zorder=3))   # the outer wall

    # The lamina: the fluid between the upper wall and the dividing streamline, which all leaves through the gap.
    div = dividing_streamline()
    lam = np.vstack([div, [[G1, 1.0 + WALL], [G0, 1.0 + WALL], [G0, 1.0], [XL, 1.0]]])
    ax.add_patch(Polygon(lam, closed=True, color=C_LAM, alpha=0.30, lw=0, zorder=1))
    ax.plot(div[:, 0], div[:, 1], color=C_INK, lw=1.1, ls=(0, (4, 2.5)), zorder=4)
    for k in (0.3, 0.7):
        p = lamina_streamline(k)
        ax.plot(p[:, 0], p[:, 1], color=C_OUT, lw=1.1, zorder=4)
        arrow_on(ax, p, 0.93, C_OUT)
    for y in (1.0 + WALL + 0.29, 1.0 + WALL + 0.37):                  # filtrate from the gaps upstream
        ax.plot([XL, XR - 0.15], [y, y], color=C_MUTED, lw=0.9, zorder=2)
        ax.add_patch(FancyArrowPatch((XR - 0.6, y), (XR - 0.15, y), arrowstyle='-|>', mutation_scale=10,
                                     color=C_MUTED, lw=0.9))
    for y in (0.30, 0.10):                                             # the center channel's flow
        p = center_streamline(y)
        ax.plot(p[:, 0], p[:, 1], color=C_STAY, lw=1.0, zorder=2)
        arrow_on(ax, p, 0.97, C_STAY)
    ax.plot([XL, XR], [0.0, 0.0], color=C_MUTED, lw=0.8, ls=(0, (8, 3, 1.5, 3)), zorder=2)

    # The two particles and the paths of their centers. The small one's center lies inside the lamina, so it follows
    # a lamina streamline through the gap. The large one's center cannot come closer to the wall than its radius,
    # which places it outside the lamina, so it passes the gap and stays in the center channel.
    k_small = A_SMALL / DELTA
    p = lamina_streamline(k_small)
    p = p[p[:, 0] >= X_SMALL]
    ax.plot(p[:, 0], p[:, 1], color=C_OUT, lw=1.3, ls=(0, (2, 2)), zorder=5)
    ax.add_patch(Circle((X_SMALL, 1.0 - A_SMALL), A_SMALL, color=C_OUT, zorder=6))
    end = p[-12]
    ax.add_patch(Circle(end, A_SMALL, color=C_OUT, alpha=0.45, zorder=6))
    yl = 1.0 - A_LARGE
    ax.plot([X_LARGE, XR - 0.5], [yl, yl], color=C_STAY, lw=1.3, ls=(0, (2, 2)), zorder=5)
    ax.add_patch(Circle((X_LARGE, yl), A_LARGE, color=C_STAY, zorder=6))
    ax.add_patch(Circle((XR - 0.75, yl), A_LARGE, color=C_STAY, alpha=0.40, zorder=6))
    for x, a, c in ((X_SMALL, A_SMALL, 'white'), (X_LARGE, A_LARGE, 'white')):
        ax.plot([x], [1.0 - a], 'o', ms=2.2, color=c, zorder=7)
    ax.annotate('', xy=(X_LARGE + 0.02, 1.0), xytext=(X_LARGE + 0.02, yl),
                arrowprops=dict(arrowstyle='<->', color='white', lw=0.9, shrinkA=0, shrinkB=0), zorder=8)
    ax.text(X_LARGE + 0.08, 1.0 - 0.55 * A_LARGE, '$a$', color='white', fontsize=10, zorder=8, va='center')

    # The lamina's width, marked where the lamina is undisturbed.
    xd = XL + 0.35
    ax.annotate('', xy=(xd, 1.0), xytext=(xd, 1.0 - DELTA),
                arrowprops=dict(arrowstyle='<->', color=C_INK, lw=0.9, shrinkA=0, shrinkB=0), zorder=8)
    ax.text(xd + 0.08, 1.0 - 0.5 * DELTA, r'$\delta$', fontsize=12, va='center', color=C_INK, zorder=8)

    # Labels sit on the walls or outside the fluid.
    ax.text(-1.6, 1.0 + 0.5 * WALL, 'post', color='white', ha='center', va='center', fontsize=9, zorder=8)
    ax.text(1.85, 1.0 + 0.5 * WALL, 'post', color='white', ha='center', va='center', fontsize=9, zorder=8)
    ax.annotate('gap', xy=(0.5 * (G0 + G1), 1.0 + WALL), xytext=(0.5 * (G0 + G1) + 0.75, 1.0 + WALL + W_SIDE + 0.27),
                ha='left', va='center', fontsize=10, color=C_INK, zorder=9,
                arrowprops=dict(arrowstyle='-', color=C_INK, lw=0.7, shrinkA=2, shrinkB=0))
    ax.text(XR + 0.1, 1.0 + WALL + 0.5 * W_SIDE, 'side channel\n(filtrate)', va='center', fontsize=10)
    ax.text(XR + 0.1, 0.52, 'center channel\n(retentate)', va='center', fontsize=10)
    ax.text(XR + 0.1, 0.0, 'center line', va='center', fontsize=9, color=C_MUTED)
    ax.text(XL, -0.13, 'flow', fontsize=9, color=C_MUTED, va='top')
    ax.add_patch(FancyArrowPatch((XL + 0.45, -0.19), (XL + 1.25, -0.19), arrowstyle='-|>', mutation_scale=10,
                                 color=C_MUTED, lw=0.9))
    ax.set_xlim(XL, XR + 1.45)
    ax.set_ylim(-0.32, 1.0 + WALL + W_SIDE + 0.42)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('(a) Plan view of one gap in the upper wall (not to scale)', fontsize=11, loc='left')


def draw_profile(ax):
    s = np.linspace(0.0, 1.0, 400)
    u = 6 * s * (1 - s)                                             # plane Poiseuille flow, mean velocity 1
    ax.plot(u, s, color=C_INK, lw=1.4)
    m = s <= X_LAM
    ax.fill_betweenx(s[m], 0, u[m], color=C_LAM, alpha=0.55, lw=0)
    ax.axhline(X_LAM, color=C_INK, lw=0.9, ls=(0, (4, 2.5)))
    ax.set_ylim(1.0, 0.0)                                           # the upper wall at the top, as in panel (a)
    ax.set_xlim(0.0, 1.75)
    ax.set_xlabel(r'velocity $u/\bar{u}$')
    ax.set_ylabel(r'distance from the upper wall, $y/w_c$')
    ax.text(0.62, X_LAM - 0.012, r'$x = \delta/w_c$', fontsize=10, va='bottom')
    ax.annotate('$f$', xy=(0.12, 0.045), xytext=(0.62, 0.20), fontsize=11, color=C_OUT, va='center',
                arrowprops=dict(arrowstyle='-', color=C_OUT, lw=0.8, shrinkA=1, shrinkB=0))
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.set_title('(b) Velocity across the center channel', fontsize=11, loc='left')


def main(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    fig = plt.figure(figsize=(12.0, 7.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 2.6], height_ratios=[1.0, 0.85], wspace=0.10, hspace=0.16)
    ax_a = fig.add_subplot(gs[0, :])
    ax_a.set_anchor('SW')
    draw_plan(ax_a)
    draw_profile(fig.add_subplot(gs[1, 0]))
    key = fig.add_subplot(gs[1, 1])
    key.axis('off')
    rows = [(C_LAM, 0.55, 's', 'the lamina, of width $\\delta$, is the layer next to the wall that leaves through the gap;\n'
                               'its share of the center flow is the shaded area in (b), $f = 3x^2 - 2x^3$ with\n'
                               '$x = \\delta/w_c$ (App. A, Result 1)'),
            (C_OUT, 1.0, 'o', 'a particle of radius $a < \\delta$ has its center inside the lamina and leaves with '
                              'the filtrate'),
            (C_STAY, 1.0, 'o', 'a particle of radius $a > \\delta$ keeps its center outside the lamina and stays in\n'
                               'the center channel, so the steric cutoff diameter is $D = 2\\delta$')]
    for j, (c, al, mk, text) in enumerate(rows):
        y = (0.86, 0.60, 0.40)[j]
        key.plot([0.06], [y], mk, ms=10 if mk == 'o' else 12, color=c, alpha=al, transform=key.transAxes)
        key.text(0.10, y, text, transform=key.transAxes, va='center', fontsize=10, linespacing=1.25)
    key.text(0.10, 0.12, 'The sketch exaggerates the lamina about eight times ($x$ = %.2f and $f$ = %.3f here). In a '
             'device,\n$x$ is about 0.01, so $f$ is a few parts in ten thousand.' % (X_LAM, F_GAP),
             transform=key.transAxes, va='center', fontsize=9.5, color=C_MUTED, linespacing=1.25)
    fig.suptitle('Skimming at one gap of a CIF device', fontsize=13, x=0.02, ha='left', y=0.985)
    path = os.path.join(out_dir, 'lamina_sketch.png')
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return path


if __name__ == '__main__':
    print(main(sys.argv[1] if len(sys.argv) > 1 else os.path.join('docs', 'figures')))
