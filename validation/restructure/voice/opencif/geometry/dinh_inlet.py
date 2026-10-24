"""opencif.geometry.dinh_inlet: the inlet of Dinh 2024 Fig. 2A drawn as depicted.

The inlet consists of the access port, the neck, a tube holding a two-column chevron filter that ends in an apex
post, a bulb crossed by a row of nine dashes (the second filter stage), and a funnel into the stem of the center
channel. Every number below was measured on the published figure (r4/meas_dinh_filter.py, r4/meas_dinh_env.py). The
inset of Fig. 2A is a native raster of 9.89 um/px, and two independent scales agree with it, since the inset spans
the 3.2 mm chip (9.85 um/px) and the head's row pitch (11.1 px) equals 2 (GS + wall) = 110 um (9.91 um/px). The
main layout of Fig. 2A, at 32.6 um/px near the inlet (bulb 36 px = 1168 um, tube 13 px = 430 um), gives the port and
the neck, which lie outside the inset.

The raster draws black (fluid) features about 0.6 px narrower and white (solid) ones 0.6 px wider than they are, so
the half-widths of fluid are corrected by +3 um and the post sizes by -6 um. The coordinate x runs along the flow
from the row-1 mouth (gap 1 at x = 0), and y is lateral. Dinh prints none of these dimensions. The filter does not
enter the CIF closure, because it is upstream of gap 1 and the closure is flow-driven, so it only has to respect the
feature floor."""
import json
import os
import numpy as np
from shapely.geometry import Point, Polygon, box
from shapely.ops import unary_union
from shapely import affinity

HERE = os.path.dirname(os.path.abspath(__file__))
_MEAS = os.path.join(HERE, '..', 'data', 'r4', 'dinh_filter_meas.json')

PORT_X, PORT_R = -5118.0, 1020.0        # main layout: port center 157 px above row 1, diameter 62.5 px
# the neck and flare (inset rows 193 to 229 and the main layout), as the half-width (um) at x (um) with the +3 um
# raster correction
NECK = [(-4108.0, 132.0), (-4020.0, 132.0), (-3870.0, 107.0), (-3822.0, 123.0), (-3793.0, 148.0),
        (-3763.0, 171.0), (-3733.0, 188.0), (-3694.0, 193.0), (-3600.0, 203.0)]
CHEV = dict(n=33, x0=-3736.0, pitch=84.2, y_top=131.0, y_bot=69.0, x_conv=-1798.0, L=79.0, W=37.0)
APEX = dict(x=-993.0, base=82.0, height=56.0)
CROSS = dict(x=-830.0, ys=(0.0, 125.0, -125.0, 250.0, -250.0, 378.0, -378.0, 505.0, -505.0), L=76.0, W=37.0)
BULB_HALF = 587.0


def _envelope():
    """Return the measured envelope, as (x, half-width) rows in um, from data/r4/dinh_filter_meas.json."""
    d = json.load(open(_MEAS))
    return np.array(d['envelope_used'], float)


def stadium(xc, yc, L, W, ang_deg, n=24):
    """Return a stadium of end-to-end length L and width W centered at (xc, yc), its long axis at ang_deg from +x."""
    r = W / 2.0
    core = box(-(L / 2.0 - r), -1e-9, L / 2.0 - r, 1e-9) if L > W else Point(0, 0)
    s = core.buffer(r, quad_segs=n)
    s = affinity.rotate(s, ang_deg, origin=(0, 0))
    return affinity.translate(s, xc, yc)


def outline(wc0):
    """Return the fluid outline of the inlet from the port to x = 1 um, where the stem continues as the trunk."""
    env = _envelope()
    xs = list(env[:, 0])
    hs = list(env[:, 1])
    # the neck, from inside the port disc (x = PORT_X + 0.9 PORT_R) through the measured hourglass and flare
    nx = [PORT_X + 0.9 * PORT_R] + [x for x, h in NECK]
    nh = [NECK[0][1]] + [h for x, h in NECK]
    # The stem blends the measured 36.3 um half-width into wc0 / 2 between x = -200 and -50 um.
    X = np.array(nx + xs + [1.0])
    H = np.array(nh + hs + [wc0 / 2.0])
    j = X > -200.0
    H[j] = H[j] + (wc0 / 2.0 - H[j]) * np.clip((X[j] + 200.0) / 150.0, 0, 1)
    up = list(zip(X, H))
    poly = Polygon([(x, -h) for x, h in up] + [(x, h) for x, h in up[::-1]]).buffer(0)
    return unary_union([poly, Point(PORT_X, 0.0).buffer(PORT_R, quad_segs=64)])


def posts():
    """Return the filter's posts, which are the two chevron columns, the apex post, and the row of dashes across the
    bulb."""
    out = []
    c = CHEV
    for k in range(c['n']):
        x = c['x0'] + k * c['pitch']
        t = min(max((x - c['x0']) / (c['x_conv'] - c['x0']), 0.0), 1.0)
        y = c['y_top'] + (c['y_bot'] - c['y_top']) * t
        out.append(stadium(x, +y, c['L'], c['W'], +45.0))     # long axis outward-downstream
        out.append(stadium(x, -y, c['L'], c['W'], -45.0))
    a = APEX
    tri = Polygon([(a['x'] - a['height'] / 2.0, 0.0), (a['x'] + a['height'] / 2.0, a['base'] / 2.0),
                   (a['x'] + a['height'] / 2.0, -a['base'] / 2.0)])
    out.append(tri.buffer(-10.0).buffer(10.0, quad_segs=16))       # rounded triangle (r 10 >= DRC d/2), points upstream
    for y in CROSS['ys']:
        out.append(stadium(CROSS['x'], y, CROSS['L'], CROSS['W'], 90.0))   # dashes, long across the bulb
    return out


def inlet(spec):
    """Return (polygon, info) of the depicted Dinh inlet in leg-1 local coordinates."""
    fl = outline(spec.wc0)
    ps = posts()
    g = fl.difference(unary_union(ps))
    info = dict(x_port_c=PORT_X, port_r=PORT_R, n_posts=len(ps), x_f0=CHEV['x0'] - 60.0, x_f1=CROSS['x'] + 60.0,
                Hf=BULB_HALF, kind='dinh', bulb_half=BULB_HALF)
    return g, info
