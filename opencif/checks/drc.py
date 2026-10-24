"""drc.py: a uniform fabrication audit of a DXF photomask, vector-based and tile-based, with a confirmation pass for
the flagged pieces.

The fluid is the union of the closed LWPOLYLINEs on every layer whose name starts with CHANNELS, which is the official
convention (the 7-7, 7-21, 7-18, and official writers all use it, and the student writer uses CHANNELS too). The
solid is the polygons on the POSTS, PILLARS, and ISLANDS layers (informational solid layers) and, for the wall test,
the band within 3d of the fluid. The tests run at the feature size d, whose default of 18.9 um catches anything
strictly narrower than the 19 um floor, and with r = d/2 they are

  open   fluid minus its morphological opening by a disk of diameter d   -> sub-d open features
  solid  band minus its opening, band = fluid.buffer(3d) - fluid          -> sub-d walls, posts, and slivers

The opening also shaves every convex corner, and those pieces (thickness ~0.34 r, extent ~1.4 r) are filtered out.
A lost region is reported only if it is a strip of near-constant width (fill >= 0.7, long side L >= 1.5 t, and
t >= 0.5 um, where t is the diameter of its maximum inscribed circle) or if t >= d/2 (see keep()). The other outputs
are the invalid (self-intersecting) rings, the overlap area among the fluid polygons (the sum of their areas minus the
area of their union; overlaps are by design in most writers, but they matter for an even-odd fill), the overlap of the
fluid with the solid layers, and the connected components of the fluid. The usage is

    python3 drc.py OUT.json file.dxf [d_um] [--islands]
"""
import sys, json, time, math, os
import ezdxf
from ezdxf import path as ezpath
import shapely
from shapely.geometry import Polygon, box
from shapely.ops import unary_union
from shapely.strtree import STRtree


def load(fname):
    """Read the closed LWPOLYLINEs of a DXF file and return the fluid polygons, the solid polygons, the number of
    invalid rings, and the number of polylines with bulges. A polyline with bulges is flattened to within 0.02 um, and
    an invalid ring is counted and then repaired with make_valid."""
    doc = ezdxf.readfile(fname)
    fluid, solid, invalid, nbulge = [], [], 0, 0
    for e in doc.modelspace().query('LWPOLYLINE'):
        lay = e.dxf.layer.upper()
        if not e.closed:
            continue
        if any(abs(b[0]) > 0 for b in e.get_points('b')):
            nbulge += 1
            pts = [(v.x, v.y) for v in ezpath.make_path(e).flattening(0.02)]
        else:
            pts = [(x, y) for x, y in e.get_points('xy')]
        if len(pts) < 3:
            continue
        pg = Polygon(pts)
        if not pg.is_valid:
            invalid += 1
            pg = shapely.make_valid(pg)
            pg = unary_union([g for g in getattr(pg, 'geoms', [pg]) if g.geom_type in ('Polygon', 'MultiPolygon')])
        if pg.is_empty:
            continue
        if lay.startswith('CHANNELS'):
            fluid.append(pg)
        elif lay in ('POSTS', 'PILLARS', 'ISLANDS'):
            solid.append(pg)
    return fluid, solid, invalid, nbulge


def pieces(g):
    """Return the parts of a geometry as a list, which holds the geometry itself when it has no parts and is empty
    for an empty geometry."""
    if g.is_empty:
        return []
    return list(g.geoms) if hasattr(g, 'geoms') else [g]


def keep(region, r, d, core):
    """Classify the lost pieces of one tile. A 'strip' is a genuine sub-d feature, either of near-constant width
    (fill >= 0.7, L >= 1.5 t, and t >= 0.5 um) or of width t >= d/2, and it is listed, whereas every other piece is a
    corner or cusp artifact of the opening, which is counted rather than listed."""
    out, ncorner = [], 0
    for p in pieces(region):
        if p.geom_type != 'Polygon' or p.area < 0.5:
            continue
        if not core.contains(p.representative_point()):
            continue                      # classified on the uncut piece and counted in the tile that holds it
        mic = shapely.maximum_inscribed_circle(p, tolerance=0.02)
        t = 2.0 * shapely.length(mic)
        rect = p.minimum_rotated_rectangle
        xs, ys = rect.exterior.coords.xy
        sides = sorted(math.hypot(xs[i + 1] - xs[i], ys[i + 1] - ys[i]) for i in range(2))
        L = sides[1]
        fill = p.area / max(L * t, 1e-9)
        if (fill >= 0.7 and L >= 1.5 * t and t >= 0.5) or t >= 0.5 * d:
            c = p.representative_point()
            out.append(dict(x=round(c.x, 1), y=round(c.y, 1), width=round(float(t), 2),
                            length=round(float(L), 1), area=round(float(p.area), 1), fill=round(float(fill), 2)))
        else:
            ncorner += 1
    return out, ncorner


def audit(fname, d=18.9, T=1500.0, subtract_islands=False):
    """Audit one DXF file at the feature size d (um) on square tiles of side T (um) and return the result as a dict.
    With subtract_islands, the solid polygons are first subtracted from the fluid."""
    t0 = time.time()
    fluid, solid, invalid, nbulge = load(fname)
    if subtract_islands:            # a writer whose ISLANDS layer holds the solid islands inside CHANNELS outlines
        isl = unary_union(solid)
        fl2 = []
        for g in fluid:
            d_ = g.difference(isl)
            fl2 += [q for q in (d_.geoms if hasattr(d_, 'geoms') else [d_]) if q.geom_type == 'Polygon' and not q.is_empty]
        fluid = fl2
    r = d / 2.0
    M = 3.2 * d
    tree = STRtree(fluid)
    stree = STRtree(solid) if solid else None
    allb = unary_union([box(*g.bounds) for g in fluid]).bounds
    x0, y0, x1, y1 = allb
    nx, ny = max(1, math.ceil((x1 - x0) / T)), max(1, math.ceil((y1 - y0) / T))
    lost_open, lost_solid = [], []
    ncorner_open = ncorner_solid = 0
    sum_area = sum(g.area for g in fluid)
    union_area = 0.0
    post_overlap = 0.0
    for ix in range(nx):
        for iy in range(ny):
            core = box(x0 + ix * T, y0 + iy * T, x0 + (ix + 1) * T, y0 + (iy + 1) * T)
            big = core.buffer(M, join_style=2)
            idx = tree.query(big)
            if len(idx) == 0:
                continue
            fl = unary_union([fluid[i] for i in idx]).intersection(big)
            union_area += fl.intersection(core).area
            if stree is not None:
                sidx = stree.query(core)
                if len(sidx):
                    sl = unary_union([solid[i] for i in sidx]).intersection(core)
                    post_overlap += sl.intersection(fl).area
            op = fl.buffer(-r, quad_segs=16).buffer(r, quad_segs=16)
            lo = fl.difference(op)
            k, c = keep(lo, r, d, core); lost_open += k; ncorner_open += c
            band = fl.buffer(3 * d, quad_segs=16).difference(fl).intersection(big)
            opb = band.buffer(-r, quad_segs=16).buffer(r, quad_segs=16)
            ls = band.difference(opb)
            # The band's own outer rim is discarded, that is, only the lost solid within d of the fluid on both sides is
            # kept.
            k, c = keep(ls, r, d, core); lost_solid += k; ncorner_solid += c
    # The confirmation pass (audit v2) is needed because a lost piece that reaches the edge of its tile's window can be
    # a clipping artifact, in which the window cut leaves a sliver of a channel that crosses it and a long sliver can
    # reach the core. The pass therefore re-opens every flagged piece on a larger window centered on it and keeps the
    # piece only if it is still lost there, away from that window's own edge. Audit v1 flagged such slivers on the
    # 3-leg Dinh devices only, whose tile grid differs from that of their linear twins with identical local geometry.
    def confirm(lst, solid_phase):
        kept, dropped = [], 0
        for q in lst:
            pad = max(6 * d, 2.0 * q['length'] + 4 * d)
            win = box(q['x'] - pad, q['y'] - pad, q['x'] + pad, q['y'] + pad)
            idx = tree.query(win.buffer(3 * d))
            fl = unary_union([fluid[i] for i in idx]).intersection(win)
            if solid_phase:
                band = fl.buffer(3 * d, quad_segs=16).difference(fl).intersection(win)
                reg = band.difference(band.buffer(-r, quad_segs=16).buffer(r, quad_segs=16))
            else:
                reg = fl.difference(fl.buffer(-r, quad_segs=16).buffer(r, quad_segs=16))
            P = shapely.geometry.Point(q['x'], q['y'])
            reg = reg.intersection(win.buffer(-3 * d))          # drop this window's own edge artifacts
            hit = [p for p in pieces(reg) if p.geom_type == 'Polygon' and p.area >= 0.5 and p.distance(P) < 1.0]
            if hit:
                kept.append(q)
            else:
                dropped += 1
        return kept, dropped
    lost_open, drop_o = confirm(lost_open, False)
    lost_solid, drop_s = confirm(lost_solid, True)
    # the connected components of the whole fluid, found cheaply from the union of all its polygons
    U = unary_union(fluid)
    comps = pieces(U)
    comp_areas = sorted((c.area for c in comps), reverse=True)
    res = dict(file=os.path.basename(fname), d_um=d, n_fluid=len(fluid), n_solid=len(solid), invalid_rings=invalid,
               bulge_polylines=nbulge, bounds=[round(v, 1) for v in allb],
               fluid_area=round(U.area, 1), overlap_area=round(sum_area - U.area, 1),
               solid_layer_overlap_with_fluid=round(post_overlap, 2),
               n_components=len(comps), largest_component_frac=round(comp_areas[0] / U.area, 6) if comps else 0,
               small_components=[round(a, 1) for a in comp_areas[1:30]],
               n_lost_open=len(lost_open), n_lost_solid=len(lost_solid), tiling_artefacts_dropped=[drop_o, drop_s], corner_artefacts_open=ncorner_open, corner_artefacts_solid=ncorner_solid,
               lost_open_min_width=min([q['width'] for q in lost_open], default=None),
               lost_solid_min_width=min([q['width'] for q in lost_solid], default=None),
               lost_open=sorted(lost_open, key=lambda q: q['width'])[:60],
               lost_solid=sorted(lost_solid, key=lambda q: q['width'])[:60],
               seconds=round(time.time() - t0, 1))
    return res


if __name__ == '__main__':
    out, f = sys.argv[1], sys.argv[2]
    d = float(sys.argv[3]) if len(sys.argv) > 3 and not sys.argv[3].startswith('--') else 18.9
    res = audit(f, d, subtract_islands=('--islands' in sys.argv))
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ('lost_open', 'lost_solid')}))
    for k in ('lost_open', 'lost_solid'):
        for q in res[k][:8]:
            print('  ', k, q)
