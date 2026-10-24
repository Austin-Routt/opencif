"""opencif.figures.drawing: draw a window of a device polygon in a chosen figure style, for side-by-side comparison
with the papers' rasters.

Coordinates are in um, and `window` = (x0, y0, x1, y1) is given in the device frame. The rotation `rot` (in degrees,
counterclockwise) is applied about `origin` before windowing, so a window can be taken in the paper's orientation.
The four styles are these:

    layout   fluid black on white (Dinh 2024 Fig. 2A CAD style)
    schem    fluid gray, solid black (Xia 2016 Fig. 1b-e schematic style)
    micro    fluid light gray, solid mid gray, dark 1-px edges (bright-field micrograph style)
    lines    white fill, black outlines of the fluid (Xia 2016 Figs. 1a and 3a, Gifford Fig. 4 schematic style)"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MPath
from shapely import affinity
from shapely.geometry import box

STYLES = dict(layout=dict(bg='white', fc='black', ec='black', lw=0.0),
              schem=dict(bg='black', fc='#b4b4b4', ec='#b4b4b4', lw=0.0),
              micro=dict(bg='#8e8e8e', fc='#cfcfcf', ec='#2a2a2a', lw=0.9),
              lines=dict(bg='white', fc='white', ec='black', lw=0.7))


def _patch(g, **kw):
    """Return one PathPatch of the polygons in g, holes included, or None when g has no polygon to draw."""
    geoms = list(g.geoms) if hasattr(g, 'geoms') else [g]
    verts, codes = [], []
    for p in geoms:
        if p.is_empty or p.geom_type != 'Polygon':
            continue
        for ring in [p.exterior] + list(p.interiors):
            xy = np.asarray(ring.coords)
            if len(xy) < 3:
                continue
            verts.extend(xy.tolist())
            codes.extend([MPath.MOVETO] + [MPath.LINETO] * (len(xy) - 2) + [MPath.CLOSEPOLY])
    if not verts:
        return None
    return PathPatch(MPath(verts, codes), **kw)


def clip(fluid, window, rot=0.0, origin=(0.0, 0.0), flipx=False):
    """Rotate the fluid by rot degrees about origin, mirror it in x when flipx is set, and cut it to the window (when
    one is given) with a 1 um margin."""
    g = fluid
    if rot:
        g = affinity.rotate(g, rot, origin=origin)
    if flipx:
        g = affinity.scale(g, -1.0, 1.0, origin=origin)
    if window is not None:
        x0, y0, x1, y1 = window
        g = g.intersection(box(x0 - 1, y0 - 1, x1 + 1, y1 + 1))
    return g


def draw(ax, fluid, window, style='layout', rot=0.0, origin=(0.0, 0.0), flipx=False, pre_clipped=False):
    """Draw the fluid in the window on ax in the given style, with equal axes set to the window, and return the clipped
    geometry."""
    st = STYLES[style]
    g = fluid if pre_clipped else clip(fluid, window, rot, origin, flipx)
    ax.set_facecolor(st['bg'])
    p = _patch(g, facecolor=st['fc'], edgecolor=st['ec'], linewidth=st['lw'])
    if p is not None:
        ax.add_patch(p)
    x0, y0, x1, y1 = window
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect('equal')
    return g


def scalebar(ax, length_um, label=None, loc=(0.06, 0.06), color='black', lw=3, fs=8):
    """Draw a horizontal scale bar of length_um at the fractional position loc of the axes, labeled with its length in
    um unless a label is given."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xa = x0 + loc[0] * (x1 - x0)
    ya = y0 + loc[1] * (y1 - y0)
    ax.plot([xa, xa + length_um], [ya, ya], color=color, lw=lw, solid_capstyle='butt')
    ax.text(xa + length_um / 2.0, ya + 0.02 * (y1 - y0), label or ('%g um' % length_um), ha='center', va='bottom',
            color=color, fontsize=fs)


def paper(ax, img, extent, cmap='gray'):
    """Show a paper's raster with its extent in um, given as (x0, x1, y0, y1)."""
    ax.imshow(img, extent=extent, cmap=cmap if img.ndim == 2 else None, interpolation='lanczos', aspect='equal')
