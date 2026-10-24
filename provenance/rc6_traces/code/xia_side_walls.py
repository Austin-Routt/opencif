"""Xia Fig. 3 panels f and g: where the image's side-channel outer walls lie against each mask's, on each panel's own
registration (the record's chamfer map, nothing re-fitted).  The walls run along the image's columns, so each is found
as the darkest column of the panel's column profile within +-30 px of the mask's wall mapped onto the image.  Offset
> 0: the image's wall lies inside the mask's side channel (the device's channel narrower than drawn).  Panel f is at its
own scale (1.774 µm/px, from its pitch); at the 1.798 of the other panels its widths would be 1.35 % larger (about 7 µm
on each side channel), which shrinks its offsets by about that much.  The leg's centre (which side is which) is the
mean y of the mask's posts in the panel's region (the posts flank the centre channel; recorded with the region's
mid-line, which the first form of this script took, for comparison).

usage: python3 -I xia_side_walls.py RASTERS RESULTS_XIA OPENCIF_ROOT VARIANTS_DIR OUT.json"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trace_lib as T

RAS, RES, ROOT, VAR, OUT = sys.argv[1:6]
from panels_xia import FIG3_SHA
im = T.Raster(os.path.join(RAS, 'xia5-000.png'), FIG3_SHA).im.astype(float)       # the Fig. 3 raster, checked by its sha
res = {}
for tag in ('rc6', 'dfill', 'fig3'):
    for panel, f in (('f', 'xia_results_fig3f_%s.json' % tag), ('g', 'xia_results_%s.json' % tag)):
        p = os.path.join(RES, f)
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        mp = r['map_chamfer']
        M = T.Map(mp['origin_px'], mp['ex_px_per_um'], mp['ey_px_per_um'])
        c0, r0, c1, r1 = r.get('micrograph_box', r['crop'])
        k = r['px_per_um']
        stem = r['mask'][:-4]
        dxf = os.path.join(ROOT, 'runs', stem, stem + '.dxf') if stem == 'xia2016_rs7_depicted' else os.path.join(VAR, stem, stem + '.dxf')
        mask = T.Mask(dxf)
        if mask.sha256 != r['mask_sha256']:
            raise SystemExit('%s: the DXF is not the one traced' % f)
        prof = im[r0 + 20:r1 - 20, c0:c1].mean(axis=0)
        walls, long_ = [], []
        for e in mask.edges(*r['region_um']):
            a = np.asarray(e, float)
            if e.tag != 'channel' or np.ptp(a[:, 0]) < 300 or np.ptp(a[:, 1]) > 20:
                continue
            col = float(np.mean(M(a)[:, 0]))
            if c0 + 5 < col < c1 - 5:
                long_.append((a, col))
        posts = [np.asarray(e, float) for e in mask.edges(*r['region_um']) if e.tag == 'post']
        y_mid = float(np.mean(np.vstack(posts)[:, 1]))       # the leg's centre: the posts' mean y (both sides of the centre channel)
        y_region = 0.5 * (r['region_um'][2] + r['region_um'][3])
        for a, col in long_:
            lo, hi = int(max(c0, col - 30)) - c0, int(min(c1 - 1, col + 30)) - c0
            j = int(np.argmin(prof[lo:hi])) + lo + c0
            centre_col = float(M(np.array([[np.mean(a[:, 0]), y_mid]]))[0, 0])
            inside = (j - col) * np.sign(centre_col - col)       # > 0: the image's wall lies towards the leg's centre
            walls.append(dict(side='left' if col < centre_col else 'right', mask_col_px=col, image_col_px=j,
                              image_inside_um=float(inside / k)))
        res['%s_%s' % (panel, tag)] = dict(mask=stem, um_per_px=1.0 / k, leg_centre_y_um=y_mid, region_mid_y_um=y_region,
                                           walls=sorted(walls, key=lambda w: w['mask_col_px']))
json.dump(dict(method=__doc__.split('\n\n')[0].replace('\n', ' '), panels=res), open(OUT, 'w'), indent=1)
for k_, v in res.items():
    print(k_, v['mask'], ['%s %+.1f' % (w['side'], w['image_inside_um']) for w in v['walls']])
