"""Measure Dinh 2024 Fig. 2A inset (native pdfimages raster dinh-001.png): hairpin fingers (white islands inside
the black hairpin outlines), arm goal-posts, trunk.  Output: per row reach, y, and arm tops, in native px."""
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
IM = np.asarray(Image.open('../r3/papers/dinh-001.png').convert('L')).astype(float)
X0, Y0, X1, Y1 = 20, 440, 370, 820
sub = IM[Y0:Y1, X0:X1]
for thr in (128, 160, 190):
    white = sub > thr
    lab, n = ndi.label(white)
    # enclosed = not touching the crop border and not the big background
    objs = ndi.find_objects(lab)
    out = []
    for i, sl in enumerate(objs):
        ys, xs = sl
        h, w = ys.stop - ys.start, xs.stop - xs.start
        area = int((lab[sl] == i + 1).sum())
        if ys.start == 0 or xs.start == 0 or ys.stop == sub.shape[0] or xs.stop == sub.shape[1]:
            continue
        out.append((i + 1, xs.start + X0, xs.stop + X0, ys.start + Y0, ys.stop + Y0, w, h, area))
    print('thr', thr, 'components', n, 'enclosed', len(out))
    big = [o for o in out if o[5] >= 15 or o[6] >= 15]
    for o in sorted(big, key=lambda o: (o[3], o[1])):
        print('  x %4d-%4d  y %4d-%4d  w %3d h %3d area %4d' % o[1:])

# ---- record the values the generator uses (meander3/gen/depicted.py, transcribed from the listing above) ----------
import json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'meander3', 'gen'))
from opencif.geometry import depicted as DP
SC, AX, ROW1 = DP.DINH_SC, DP.DINH_AX, DP.DINH_ROW1_PX
rec = dict(scale_um_per_px=SC, axis_px=AX, row1_mouth_px=ROW1,
           note='arm tops: outer edge of the black U-turn at the inner lane column; finger extents at threshold 160; '
                'the raster thins fluid by ~3 um per side (corrected below)',
           arm_top_px=DP.DINH_ARM_TOP_PX,
           arm_top_outer_edge_mm_above_row1=[round(-((t - ROW1) * SC - 3.0) / 1000.0, 3) for t in DP.DINH_ARM_TOP_PX],
           finger_left_px=DP.DINH_FIN_L, finger_right_px=DP.DINH_FIN_R,
           finger_reach_um=[round(0.5 * ((r - 0.5 - AX) + (AX - (l - 0.5))) * SC - 3.0, 1)
                            for l, r in zip(DP.DINH_FIN_L, DP.DINH_FIN_R)])
json.dump(rec, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dinh_inset_meas.json'), 'w'), indent=1)
print(json.dumps(rec, indent=1))
