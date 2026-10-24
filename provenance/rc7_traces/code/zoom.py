"""inspection zooms: python3 zoom.py RESULTS.json KEY RASTER DXF c0 r0 c1 r1 OUT.png [scale]"""
import json, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trace_lib as T
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
res, key, ras, dxf = sys.argv[1:5]
c0, r0, c1, r1 = map(int, sys.argv[5:9]); outp = sys.argv[9]; sc = float(sys.argv[10]) if len(sys.argv) > 10 else 6
R = json.load(open(res)); R = R[key] if key != '-' else R
mp = R['map_chamfer']; M = T.Map(mp['origin_px'], mp['ex_px_per_um'], mp['ey_px_per_um'])
im = T.Raster(ras).im
mask = T.Mask(dxf)
# window corners -> mask box
corners = M.inverse(np.array([[c0, r0], [c1, r0], [c0, r1], [c1, r1]], float))
x0, x1 = corners[:, 0].min(), corners[:, 0].max(); y0, y1 = corners[:, 1].min(), corners[:, 1].max()
segs = [M(e) for e in mask.edges(x0, x1, y0, y1, post_max=(1e9, 1e9))]
W, H = c1 - c0, r1 - r0
fig = plt.figure(figsize=(W * sc / 100, H * sc / 100), dpi=100)
ax = fig.add_axes([0, 0, 1, 1])
ax.imshow(im[r0:r1, c0:c1], cmap='gray', vmin=0, vmax=255, interpolation='nearest', extent=(c0 - .5, c1 - .5, r1 - .5, r0 - .5))
ax.add_collection(LineCollection(segs, colors='#00d5ff', linewidths=0.8))
ax.set_xlim(c0 - .5, c1 - .5); ax.set_ylim(r1 - .5, r0 - .5); ax.axis('off')
fig.savefig(outp); print(outp)
