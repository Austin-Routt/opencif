"""full_device_composite.png: the four shipped devices at ONE scale, assembled from their whole-device strips.

full_device.py draws each device in strips at 2 µm/px (the drawing area of a strip is a whole number of pixels with x0 at
its left edge and y0 at its bottom edge, below a caption band of about 70 px).  Drawing all four devices again in one figure needs
more memory than the container has (it was stopped at 6.3 GB), so this script reduces each strip's drawing area by an
integer factor k with a box average (PIL Image.reduce), after padding it on the right and at the top to a multiple of k
so that x0 and y0 stay on pixel edges: the composite is at exactly k × 2 µm/px (V10), and nothing is re-fitted or
re-drawn.  Each row carries its device, strip and x range.
usage: python3 -I composite_from_strips.py FULL_DEVICE_DIR [K, default 5] [NAMES...]"""
import json
import math
import os
import sys

from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None
D = sys.argv[1]
K = int(sys.argv[2]) if len(sys.argv) > 2 else 5
SHIP = sys.argv[3:] or ['dinh2024_v1_depicted', 'gifford2014_fig4_design', 'gifford2014_fig6_design', 'xia2016_rs7_depicted']
HEAD = 70                                    # full_device.py: the caption band above each strip's drawing area
idx = json.load(open(os.path.join(D, 'full_device_index.json')))
try:
    font = ImageFont.truetype('DejaVuSans.ttf', 28)
except OSError:
    font = ImageFont.load_default()


def strip_boxes(dev, s):
    """full_device.render_device's strip boxes, recomputed from the index (the same formula)"""
    X0, X1, Y0, Y1 = dev['bbox_um']
    pad = 0.04 * (Y1 - Y0) + 150.0
    n = dev['strips']
    seg = (X1 - X0) / n
    bx = [(X0 + i * seg - 50.0, X0 + (i + 1) * seg + 50.0, Y0 - pad, Y1 + pad) for i in range(n)]
    return [(x0, x0 + math.ceil((x1 - x0) / s) * s, y0, y0 + math.ceil((y1 - y0) / s) * s) for x0, x1, y0, y1 in bx]


rows, meta = [], []
for name in SHIP:
    dev = idx['devices'][name]
    s = float(dev['strip_um_per_px'])
    for i, (x0, x1, y0, y1) in enumerate(strip_boxes(dev, s)):
        W, H = int(round((x1 - x0) / s)), int(round((y1 - y0) / s))
        im = Image.open(os.path.join(D, '%s_strip%d.png' % (name, i + 1))).convert('RGB')
        # the saved strip can be a pixel short of W x (H + HEAD) (the figure size in inches times 100 dpi, truncated):
        # the drawing area is the bottom H rows (its axes are anchored at the figure's bottom left), the caption band the rest
        assert im.size[0] in (W, W - 1) and im.size[1] in (H + HEAD, H + HEAD - 1), (name, i, im.size, (W, H + HEAD))
        assert im.size[0] == W, (name, i, 'width', im.size[0], W)
        area = im.crop((0, im.size[1] - H, W, im.size[1]))
        Wp, Hp = K * math.ceil(W / K), K * math.ceil(H / K)
        pad = Image.new('RGB', (Wp, Hp), (255, 255, 255))
        pad.paste(area, (0, Hp - H))                          # x0 at the left edge, y0 at the bottom edge
        small = pad.reduce(K)
        rows.append(small)
        meta.append(dict(device=name, strip=i + 1, of=dev['strips'], box_um=[x0, x0 + Wp * s, y0, y0 + Hp * s],
                         um_per_px=[(Wp * s) / small.size[0], (Hp * s) / small.size[1]]))
        del im, area, pad
gap, band = 24, 44
Wc = max(r.size[0] for r in rows)
Hc = sum(r.size[1] + band + gap for r in rows) + gap
C = Image.new('RGB', (Wc, Hc), (255, 255, 255))
dr = ImageDraw.Draw(C)
y = gap
for r, m in zip(rows, meta):
    dr.text((4, y + 6), '%s — strip %d of %d, x %.1f to %.1f mm' % (m['device'], m['strip'], m['of'], m['box_um'][0] / 1000.0,
                                                                     m['box_um'][1] / 1000.0), fill=(0, 0, 0), font=font)
    y += band
    C.paste(r, (0, y))
    y += r.size[1] + gap
scale = s * K
# a 10 mm bar at the bottom left, at the composite's scale
bar = int(round(10000.0 / scale))
dr.rectangle((20, Hc - 14, 20 + bar, Hc - 6), fill=(214, 40, 40))
dr.text((20 + bar + 10, Hc - 40), '10 mm', fill=(214, 40, 40), font=font)
out = os.path.join(D, 'full_device_composite.png')
C.save(out, optimize=True)
ok = all(abs(m['um_per_px'][0] - scale) < 1e-9 and abs(m['um_per_px'][1] - scale) < 1e-9 for m in meta)
idx['composite'] = dict(um_per_px=scale, method=('the four shipped devices\' strips (%.0f µm/px) reduced by %d by box averaging, '
                                                 'padded right and top so that each strip\'s x0 and y0 stay on pixel edges '
                                                 '(code/composite_from_strips.py)' % (s, K)),
                        rows=meta, V10=dict(ok=ok, um_per_px=scale), files=['full_device_composite.png'])
json.dump(idx, open(os.path.join(D, 'full_device_index.json'), 'w'), indent=1)
print('wrote %s: %d × %d px at %.1f µm/px, %d rows, V10 %s, %.1f MB' % (out, Wc, Hc, scale, len(rows), 'ok' if ok else 'FAILED',
                                                                      os.path.getsize(out) / 1e6))
sys.exit(0 if ok else 2)
