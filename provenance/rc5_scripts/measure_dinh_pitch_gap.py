"""Dinh 2024 Fig. 2B fourth panel (raster 982 x 1940): the post pitch and the last gap, from the posts' dark outlines.
Columns 805-828 of the right post column and 705-728 of the left (both clear of the caption arrows, which cross the
posts only at rows ~1640-1650 and ~1685-1695 near their outer ends) are averaged per row; every outline row is a local
minimum of that profile, refined to sub-pixel by a parabola.  Pill k's centre = mid-point of its top and bottom
outline; the gap between the last pill and the long wall = from the last pill's bottom outline to the long wall's top
outline.  Two caption arrows cross the posts (rows ~1646 and ~1688); their shafts are dropped.  Usage: python3 measure_dinh_pitch_gap.py <raster.png>"""
import sys, json
import numpy as np
from PIL import Image
from scipy.signal import find_peaks
IM = np.asarray(Image.open(sys.argv[1]).convert('L')).astype(float)
UM = 50.0 / 46.0                                      # the panel's 50 um bar is 46 px (measure_dinh_flare.py)
out = {}
for name, (xa, xb) in (('right_column', (805, 829)), ('left_column', (705, 729))):
    prof = IM[1495:1725, xa:xb].mean(axis=1)
    pk, _ = find_peaks(-prof, prominence=15, distance=8)
    rows = []
    for i in pk:
        a, b, c = prof[i - 1], prof[i], prof[i + 1]
        den = a - 2 * b + c
        rows.append(1495 + i + (0.5 * (a - c) / den if den > 0 else 0.0))
    rows = np.array(rows)
    arrow = ((rows > 1640) & (rows < 1651)) | ((rows > 1684) & (rows < 1693))   # the caption arrows' shafts (see the
    out.setdefault('arrow_rows_dropped', {})[name] = [round(float(r), 2) for r in rows[arrow]]  # panel), not outlines
    rows = rows[~arrow]
    # outlines alternate top/bottom; the last one is the long wall's top
    tops, bots = rows[0:-1:2], rows[1:-1:2]
    n = min(len(tops), len(bots))
    cen = 0.5 * (tops[:n] + bots[:n])
    pitch_px = np.diff(cen)
    out[name] = dict(outline_rows=[round(float(r), 2) for r in rows], pill_centres=[round(float(c), 2) for c in cen],
                     pitch_px=[round(float(p), 2) for p in pitch_px], pitch_um_mean=round(float(np.mean(pitch_px) * UM), 2),
                     post_along_flow_um=round(float(np.mean(bots[:n] - tops[:n]) * UM), 2),
                     last_gap_um=round(float((rows[-1] - bots[n - 1]) * UM), 2),
                     last_gap_centre_row=round(float(0.5 * (rows[-1] + bots[n - 1])), 2),
                     last_pill_bottom_row=round(float(bots[n - 1]), 2), long_wall_top_row=round(float(rows[-1]), 2))
json.dump(out, open('dinh_pitch_gap_measurement.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
