"""meas_gifford_inset.py -- Gifford 2014 Fig. 4's bright-field inset (giff-004.png, frame x 197-829, y 52-552; bar
300 um = 107 px, 2.804 um/px): the widths across leg 2 and the end of leg 4, and the wall between them.

Method: rotate the crop by -45 deg (the angle that makes the column profile sharpest), so the legs run vertically;
average three 200-row windows column by column; the channel walls and the post outlines are dark lines, i.e. dips of
the profile.  Each feature is the strongest dip within +-6 px of its approximate position (read once from the rotated
image, r4/g4_inset_rot.png).  Widths are dip-to-dip (line centre to line centre), +-2 px (+-1 px per line).
Check on a known width: leg 4's side channels at its outlet end must read the design's 300 um."""
import json
import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.signal import find_peaks

SC = 300.0 / 107.0
IMG = np.asarray(Image.open('../r3/papers/giff-004.png').convert('L')).astype(float)
crop = IMG[52:552, 197:829]
ANG = -45.0
r = ndimage.rotate(crop, ANG, reshape=True, order=1, cval=np.nan)
APPROX = dict(leg2_wall_out=153, leg2_post_out_L=197, leg2_post_in_L=217, leg2_post_in_R=253, leg2_post_out_R=274,
              leg2_wall_in=319, leg4_wall_in=382, leg4_post_out_L=490, leg4_post_out_R=568, leg4_wall_out=674)
pos = {k: [] for k in APPROX}
for (y0, y1) in ((250, 450), (350, 550), (420, 620)):
    ok = np.sum(np.isfinite(r[y0:y1, :]), axis=0) == (y1 - y0)
    prof = np.where(ok, np.nanmean(np.where(np.isfinite(r[y0:y1, :]), r[y0:y1, :], 0.0), axis=0)
                    * (y1 - y0) / np.maximum(np.sum(np.isfinite(r[y0:y1, :]), axis=0), 1), np.nan)
    base = np.nanmedian(prof)
    pk, pr = find_peaks(-np.nan_to_num(prof, nan=base), prominence=2)
    for k, x0 in APPROX.items():
        c = [(pr['prominences'][i], int(p)) for i, p in enumerate(pk) if abs(p - x0) <= 6]
        if c:
            pos[k].append(max(c)[1])
P = {k: float(np.median(v)) for k, v in pos.items() if v}

# row pitch: autocorrelation of each post row's column profile along the leg; the 10-period peak, divided by 10
LAM = []
for (c0, c1, y0, y1) in ((492, 508, 380, 700), (548, 564, 380, 700), (199, 215, 150, 600), (255, 271, 150, 600)):
    seg = r[y0:y1, c0:c1]
    p = np.nanmean(seg, axis=1)
    p = p[np.isfinite(p)] - np.nanmean(p)
    ac = np.correlate(p, p, 'full')[len(p) - 1:]
    k1 = 10 + int(np.argmax(ac[10:24]))                  # one period (the 2-period peak is not taken)
    k10 = 10 * k1 - 8 + int(np.argmax(ac[10 * k1 - 8:10 * k1 + 9]))
    y_, y0_, y2_ = ac[k10], ac[k10 - 1], ac[k10 + 1]
    LAM.append((k10 + 0.5 * (y0_ - y2_) / (y0_ - 2 * y_ + y2_)) / 10.0 * SC)


def w(a, b):
    return round((P[b] - P[a]) * SC, 1)


out = dict(scale_um_per_px=SC, angle_deg=ANG, positions_px=P, n_windows={k: len(v) for k, v in pos.items()},
           uncertainty_um=round(2 * SC, 1),
           leg2_side_L_um=w('leg2_wall_out', 'leg2_post_out_L'), leg2_side_R_um=w('leg2_post_out_R', 'leg2_wall_in'),
           leg2_centre_um=w('leg2_post_in_L', 'leg2_post_in_R'),
           wall_between_legs_um=w('leg2_wall_in', 'leg4_wall_in'),
           leg4_side_L_um=w('leg4_wall_in', 'leg4_post_out_L'), leg4_side_R_um=w('leg4_post_out_R', 'leg4_wall_out'),
           row_pitch_um=[round(v, 2) for v in LAM], row_pitch_um_mean=round(float(np.mean(LAM)), 2),
           post_across_flow_um_line_centres=[w('leg2_post_out_L', 'leg2_post_in_L'), w('leg2_post_in_R', 'leg2_post_out_R')],
           check='leg 4 side channels at the outlet end: design 300 um')
json.dump(out, open('gifford_inset_meas.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
