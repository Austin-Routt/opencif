"""meas_xia_fig3g.py -- Xia 2016 Fig. 3g (the device's outlet end, bright field; r3/x3g.png, native raster at the
panel-i bar: 250 um = 128 px, 1.953 um/px; flow downward).  Measures:
  - the post pitch Lambda: autocorrelation of each post column's profile along the flow, 10-period peak / 10;
  - the post length across the flow, l_g, and the widths of the left side channel and the centre channel, below the
    posts (the long walls of the outlet stretch) and through the posts: distances between dark outline centres.
Writes xia_fig3g_meas.json."""
import json
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

SC = 250.0 / 128.0
IM = np.asarray(Image.open('../r3/x3g.png').convert('L')).astype(float)
out = dict(scale_um_per_px=SC)
LAM = []
for c0, c1 in ((478, 506), (600, 628)):
    p = IM[5:385, c0:c1].mean(axis=1)
    p = p - p.mean()
    ac = np.correlate(p, p, 'full')[len(p) - 1:]
    k1 = 20 + int(np.argmax(ac[20:40]))
    k10 = 10 * k1 - 10 + int(np.argmax(ac[10 * k1 - 10:10 * k1 + 11]))
    y_, y0_, y2_ = ac[k10], ac[k10 - 1], ac[k10 + 1]
    LAM.append((k10 + 0.5 * (y0_ - y2_) / (y0_ - 2 * y_ + y2_)) / 10.0 * SC)
out['post_pitch_um'] = [round(v, 2) for v in LAM]


def dips(y0, y1):
    prof = gaussian_filter1d(IM[y0:y1, :].mean(axis=0), 1.0)
    pk, pr = find_peaks(-prof, prominence=6)
    return [int(v) for v in pk], prof


# below the posts (rows 450-700): left outer wall line, then the two long walls (each a dark-outlined block)
pk, prof = dips(450, 700)
out['below_posts_dips_px'] = pk
left_wall = min(p for p in pk if 100 < p < 150)
blk = [p for p in pk if 460 < p < 650]
out['side_channel_left_um'] = round((blk[0] - left_wall) * SC, 1)
inner = [p for p in blk if 505 < p < 600]
out['centre_channel_um'] = round((max(inner) - min(inner)) * SC, 1) if len(inner) >= 2 else None
# through the posts (pill centre rows): the pill's outline centres across the flow
pk2, _ = dips(100, 360)
out['post_rows_dips_px'] = pk2
lp = [p for p in pk2 if 460 < p < 525]
rp = [p for p in pk2 if 585 < p < 650]
out['post_across_flow_um'] = [round((max(v) - min(v)) * SC, 1) for v in (lp, rp) if len(v) >= 2]
out['note'] = ('distances between dark outline centres (+-1 px per edge, +-3.9 um on a width); the right side channel '
               'runs past the crop')
json.dump(out, open('xia_fig3g_meas.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
