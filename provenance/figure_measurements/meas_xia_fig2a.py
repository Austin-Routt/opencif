"""meas_xia_fig2a.py -- Xia 2016 Fig. 2a(ii), the RS-7 schematic with its own 5 mm scale bar (page image xia-001.png,
2850 x 2130, from pdfimages of the project's copy).  A second whole-device layout of RS-7, independent of Fig. 3a: it
checks Fig. 3a's scale (27.93 um/px from its bar; the rejected reading was 31.25).  Fig. 2a(ii) is a schematic whose
outlet ring and ports are drawn smaller than Fig. 3a's, so the comparison uses the paperclip body: the top loop's outer
line to the outer line of the outermost bottom loop, measured the same way on both figures.  (The overall length to the
lowest port is recorded too, but it is not comparable.)
Panel a(ii) occupies columns ~560-1080, rows ~10-780.  The bar is the dark horizontal run at rows 589-591 near
column 880; the device's extent is that of the dark (near-black) pixels in columns 690-870 (the green callout lines
and labels are excluded by colour).  Writes xia_fig2a_meas.json."""
import json, os
import numpy as np
from PIL import Image

SP = os.environ.get('CIF_ROOT', '.')
rgb = np.asarray(Image.open(SP + '/r3/papers/xia-001.png').convert('RGB')).astype(int)
gray = rgb.mean(axis=2)
bars = []
for y in range(575, 605):
    d = np.nonzero(gray[y, 840:960] < 100)[0]
    if len(d) > 20:
        bars.append(int(d.max() - d.min() + 1))
bar_px = float(np.median(bars))
um_px = 5000.0 / bar_px
sub = rgb[10:780, 690:870]
dark = sub.sum(axis=2) < 300
rows = np.nonzero(dark.any(axis=1))[0]
cols = np.nonzero(dark.any(axis=0))[0]
length_px = float(rows.max() - rows.min() + 1)
width_px = float(cols.max() - cols.min() + 1)
out = dict(bar_px=bar_px, um_per_px=um_px, length_px=length_px, length_mm=length_px * um_px / 1000.0,
           length_rel_uncertainty=1.0 / bar_px, width_px=width_px, width_mm=width_px * um_px / 1000.0,
           rows=[int(rows.min()) + 10, int(rows.max()) + 10], cols=[int(cols.min()) + 690, int(cols.max()) + 690],
           note='overall length: the top loop\'s outer wall to the bottom of the lowest port (as Fig. 3a\'s extent); '
                'uncertainty dominated by the bar (+-1 px of 47)')
# the paperclip body on Fig. 2a(ii): the top loop's outer line (topmost dark row) to the outermost bottom loop's outer
# line (lowest dark row in columns 785-815, right of the outlet ring)
loops = np.nonzero((rgb[600:760, 785:815].sum(axis=2) < 300).any(axis=1))[0] + 600
body_px = float(loops.max() - (rows.min() + 10))
out['body_px'] = body_px
out['body_mm'] = body_px * um_px / 1000.0
# the same span on Fig. 3a (r3/x3a.png, 27.93 um/px): the top loop's outer wall (row 127, as meas_xia_fig3a.py) to the
# outer wall of the outermost bottom loop at its centre column (the last dark run before the ring, outer edge - 2 px)
IM3 = np.asarray(Image.open(SP + '/r3/x3a.png').convert('L')).astype(float)
M3 = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'meas_xia_fig3a.json')))
c5 = int(round(M3['loops']['5']['centre_x_px']))
yc = int(round(M3['loop_centre_row_bottom']['5']))
col = IM3[:, c5 - 2:c5 + 3].mean(axis=1)
dk = np.nonzero(col[yc:yc + 300] < 128)[0] + yc
runs = []
for v in dk:
    if runs and v - runs[-1][1] <= 1:
        runs[-1][1] = int(v)
    else:
        runs.append([int(v), int(v)])
gaps = [runs[i + 1][0] - runs[i][1] for i in range(len(runs) - 1)]
k_last = next(i for i, g in enumerate(gaps) if g > 40)        # the white gap between the loops and the outlet ring
bottom3 = runs[k_last][1] - 2.0
body3_px = bottom3 - 127.0
out['fig3a_body_mm'] = body3_px * M3['um_per_px'] / 1000.0
out['fig3a_body_mm_at_31_25'] = body3_px * 31.25 / 1000.0
out['body_rel_to_fig3a'] = out['body_mm'] / out['fig3a_body_mm'] - 1
out['overall_rel_to_fig3a'] = out['length_mm'] / (M3['extent_um']['length'] / 1000.0) - 1
# the outlet on Fig. 2a(ii): circle fits to the dark pixels (Kasa least squares, iterated on a +-2.5 px annulus) for the
# two ports; the ring's two wall lines lie 30-36.5 px from its centre (radial histogram), so its centreline is taken
# at their mid-radius.  Fig. 3a's values from meas_xia_fig3a.json.
ysd, xsd = np.nonzero(rgb[640:760, 700:800].sum(axis=2) < 300)
Pd = np.c_[xsd + 700.0, ysd + 640.0]


def _fit(c):
    c = np.array(c, float)
    for _ in range(8):
        sel = Pd[np.abs(np.hypot(Pd[:, 0] - c[0], Pd[:, 1] - c[1]) - c[2]) < 2.5]
        A = np.c_[2 * sel[:, 0], 2 * sel[:, 1], np.ones(len(sel))]
        q, *_ = np.linalg.lstsq(A, (sel ** 2).sum(axis=1), rcond=None)
        c = np.array([q[0], q[1], np.sqrt(q[2] + q[0] ** 2 + q[1] ** 2)])
    return c


cp2, lp2 = _fit((745, 691, 14)), _fit((745, 726, 13))
cx, cy = 745.0, 701.0
dd = np.hypot(Pd[:, 0] - cx, Pd[:, 1] - cy)
ang = np.degrees(np.arctan2(Pd[:, 1] - cy, Pd[:, 0] - cx))
side = ((np.abs(ang) < 40) | (np.abs(np.abs(ang) - 180) < 40)) & (dd > 28) & (dd < 38)   # the ring band (30-36.5 px)
ring_mid_px = 0.5 * (np.percentile(dd[side], 5) + np.percentile(dd[side], 95))
f3 = M3['um_per_px']
out['outlet'] = dict(centre_port_r_mm=cp2[2] * um_px / 1000.0, lower_port_r_mm=lp2[2] * um_px / 1000.0,
                     port_spacing_mm=abs(lp2[1] - cp2[1]) * um_px / 1000.0, ring_centreline_r_mm=ring_mid_px * um_px / 1000.0,
                     fig3a_centre_port_r_mm=M3['centre_port']['r_um'] / 1000.0,
                     fig3a_filtrate_port_r_mm=M3['filtrate_port']['r_um'] / 1000.0,
                     fig3a_port_spacing_mm=M3['ring']['filtrate_port_beyond_centre_port_um'] / 1000.0,
                     fig3a_ring_centreline_r_mm=M3['ring']['R_centreline_um'] / 1000.0)
o_ = out['outlet']
o_['rel'] = dict(centre_port_r=o_['centre_port_r_mm'] / o_['fig3a_centre_port_r_mm'] - 1,
                 filtrate_port_r=o_['lower_port_r_mm'] / o_['fig3a_filtrate_port_r_mm'] - 1,
                 port_spacing=o_['port_spacing_mm'] / o_['fig3a_port_spacing_mm'] - 1,
                 ring_centreline_r=o_['ring_centreline_r_mm'] / o_['fig3a_ring_centreline_r_mm'] - 1)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'xia_fig2a_meas.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
