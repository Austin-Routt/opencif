"""opencif.geometry.dxf: the DXF writer (R2000, in um).

Layer CHANNELS holds the fluid outlines and layer ISLANDS holds the solid islands inside them, so the fluid is
CHANNELS minus ISLANDS. Every curve is a fine polyline, with arcs drawn at no more than 3 degrees per segment. The file
carries fixed metadata (no creation time and no random GUIDs), so the same geometry gives the same bytes, and a mask
can be pinned by its sha256 (the record's `files.dxf_sha256`)."""
import hashlib
import ezdxf


def write_dxf(fluid, path, extra=None):
    """Write the fluid polygon, and any extra layers, to a DXF at path and return the number of CHANNELS and ISLANDS
    polylines. The fixed-metadata option of ezdxf is set for the write and restored afterward, so the bytes are
    reproducible without changing the option for other callers."""
    prev = ezdxf.options.write_fixed_meta_data_for_testing
    # The creation and update stamps and the GUIDs are fixed, so the bytes are reproducible.
    ezdxf.options.write_fixed_meta_data_for_testing = True
    try:
        return _write(fluid, path, extra)
    finally:
        ezdxf.options.write_fixed_meta_data_for_testing = prev


def _write(fluid, path, extra):
    doc = ezdxf.new('R2000')
    doc.header['$INSUNITS'] = 13          # micrometers
    msp = doc.modelspace()
    geoms = list(fluid.geoms) if hasattr(fluid, 'geoms') else [fluid]
    n = 0
    for g in geoms:
        msp.add_lwpolyline(list(g.exterior.coords)[:-1], close=True, dxfattribs={'layer': 'CHANNELS'})
        n += 1
        for ring in g.interiors:
            msp.add_lwpolyline(list(ring.coords)[:-1], close=True, dxfattribs={'layer': 'ISLANDS'})
            n += 1
    for lay, polys in (extra or {}).items():
        for p in polys:
            msp.add_lwpolyline(list(p.exterior.coords)[:-1], close=True, dxfattribs={'layer': lay})
    doc.saveas(path)
    return n


def sha256(path):
    """Return the sha256 of a written DXF, which is reproducible because write_dxf fixes the metadata."""
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()
