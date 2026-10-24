"""opencif, a target-first generator of Controlled Incremental Filtration (CIF) devices.

The inputs are the size cutoff D and the filtrate fraction F, which are decoupled, together with the process, the chip
area, and a choice of parts (the layout, the taper of the center channel, the meander head and its style, the inlet,
and the outlet). One call builds a device:

    import opencif
    r = opencif.generate(4.0, 0.5, chip=(12, 75), meander='pack', inlet='filter', outlet='teardrop')

The call is api.generate, and docs/OPTIONS.md lists every part with its code. The same build runs from a YAML spec
with `python -m opencif build SPEC.yaml`. The output is a DXF photomask, multipanel renders, and a validation record
that measures the drawn mask and reports, row by row, the cutoff that the mask realizes and its filtrate fraction F.

The package is organized in the parts below, and docs/ARCHITECTURE.md walks through a build:

    opencif.physics    the kernels of Appendix A (kernels), the row schedule with its invariant lamina, march, and
                       exact closure (schedule), the 1-D flow network (network), and the model and fabrication
                       scenarios (scenarios)
    opencif.geometry   the parts' shapes (shapes), the layout of legs, heads, turns, inlets, and outlets (layout), the
                       published devices' depicted parts (depicted, dinh_inlet), and the DXF writer (dxf)
    opencif.checks     the measurement of the drawn polygon (measure), the design-rule check of the DXF (drc), the
                       physics auditor (auditor), which reads the photomask back and checks P0-P10 on it
                       (docs/PHYSICS_AUDIT.md), and the design code's reports and acceptance validators (design_code)
    opencif.figures    the proof figures of each build
    opencif.spec, api, pipeline, build, compose   the specification, generate(), and the build itself
    opencif.teaching   a small, self-contained version of the whole path for a straight ladder, written to be read

The kernels are the design code's (reference/CIF_v2.py v2.1.0, which is kept byte-frozen because Appendix A names
it). opencif.physics.kernels restates them, and tests/test_kernels.py proves that they return the same numbers bit for
bit.
"""
__version__ = '1.0.0'


def generate(cutoff_um, filtrate_fraction=None, **kw):
    """Build a validated photomask, its figures, and its record from D and F (see opencif.api.generate)."""
    from .api import generate as _g
    return _g(cutoff_um, filtrate_fraction, **kw)


def build_spec(cutoff_um, filtrate_fraction=None, **kw):
    """Return the spec dict that generate() would build, without building it (see opencif.api.build_spec)."""
    from .api import build_spec as _b
    return _b(cutoff_um, filtrate_fraction, **kw)
