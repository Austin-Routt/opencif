"""Manual edits after refactor_v1.py: the package's entry module, the code stamp, paths that named the old layout.
usage: python -I patch_v1.py TREE"""
import os
import sys

T = sys.argv[1]
n_done = 0


def edit(rel, old, new, count=1):
    global n_done
    p = os.path.join(T, rel)
    s = open(p, encoding='utf-8').read()
    assert s.count(old) == count, (rel, s.count(old), old[:80])
    s = s.replace(old, new)
    open(p, 'w', encoding='utf-8').write(s)
    n_done += 1


# 1. the package entry: no sys.path changes, the version, the public functions
edit('opencif/__init__.py', '''The physics engine is the clean-room replica generator (meander3: exact gap-drop closure, exact outlet closure,
drawable sealed meander heads), with the straight-wall rule made resistance-exact; the kernels (duct law, cubic
skim kernel, kappa pairing) are imported from CIF_v2.py v2.1.0, kept byte-frozen in reference/.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.path.join(HERE, 'engine')
if ENGINE not in sys.path:
    sys.path.insert(0, ENGINE)
__version__ = '1.0.0rc8'
''', '''The package, part by part (docs/ARCHITECTURE.md walks through a build):

    opencif.physics    the kernels of Appendix A (kernels), the row schedule: invariant lamina, march and exact closure
                       (schedule), the 1-D flow network (network), the model and fabrication scenarios (scenarios)
    opencif.geometry   the parts' shapes (shapes), the layout of legs, heads, turns, inlets and outlets (layout), the
                       published devices' depicted parts (depicted, dinh_inlet) and the DXF writer (dxf)
    opencif.checks     the measurement of the drawn polygon (measure), the design-rule check of the DXF (drc), the
                       physics auditor (auditor), which reads the photomask back and checks P0-P10 on it
                       (docs/PHYSICS_AUDIT.md), and the design code's reports and acceptance validators (design_code)
    opencif.figures    the proof figures of each build
    opencif.spec, api, pipeline, build, compose   the specification, generate(), and the build itself
    opencif.teaching   a small, self-contained version of the whole path for a straight ladder, written to be read

The kernels are the design code's (reference/CIF_v2.py v2.1.0, kept byte-frozen because Appendix A names it):
opencif.physics.kernels restates them and tests/test_kernels.py proves they return the same numbers bit for bit.
"""
__version__ = '1.0.0'
''')

# 2. the DRC script's new place; the kernel pair's source
edit('opencif/build.py', "DRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'engine', 'drc_audit.py')",
     "DRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'checks', 'drc.py')")
edit('opencif/build.py', "source='CIF_v2 v2.1.0 Pairing (byte-frozen, opencif/reference); Hill fit: glossary macros')",
     "source='opencif.physics.kernels.Pairing (= CIF_v2 v2.1.0 Pairing, reference/CIF_v2.py); Hill fit: glossary macros')")

# 3. the code stamp: every Python file of the package and its data
edit('opencif/pipeline.py', '''    """sha256 over the package's source (engine included): every record says which code produced it"""
    import hashlib
    import glob
    root = os.path.dirname(os.path.abspath(__file__))
    h = hashlib.sha256()
    files = sorted(glob.glob(os.path.join(root, '*.py')) + glob.glob(os.path.join(root, 'engine', '*.py')) +
                   glob.glob(os.path.join(root, 'reference', '*.py')) + glob.glob(os.path.join(root, 'data', '*.json')) +
                   glob.glob(os.path.join(root, 'data', '*', '*.json')))
    for f in files:
        h.update(os.path.relpath(f, root).encode())
        h.update(open(f, 'rb').read())''', '''    """sha256 over the package's source: every Python file under opencif/ and its data (JSON), each hashed with its
    path relative to the package, '/' as the separator on every system.  Every record says which code produced it."""
    import hashlib
    import glob
    root = os.path.dirname(os.path.abspath(__file__))
    h = hashlib.sha256()
    files = [f for f in glob.glob(os.path.join(root, '**', '*.py'), recursive=True) +
             glob.glob(os.path.join(root, 'data', '*.json')) + glob.glob(os.path.join(root, 'data', '*', '*.json'))
             if '__pycache__' not in f]
    posix = lambda f: os.path.relpath(f, root).replace(os.sep, '/')
    files.sort(key=posix)                                # the same order on every system
    for f in files:
        h.update(posix(f).encode())
        h.update(open(f, 'rb').read())''')
# the tests' hash: the same rule ('/' and the same order on every system; the RC8 form hashed os.sep)
edit('opencif/pipeline.py', '''    tests = sorted(glob.glob(os.path.join(root, 'tests', '*.py')))
    th = hashlib.sha256()
    for f in tests:
        th.update(os.path.relpath(f, root).encode())''', '''    posix = lambda f: os.path.relpath(f, root).replace(os.sep, '/')
    tests = sorted(glob.glob(os.path.join(root, 'tests', '*.py')), key=posix)
    th = hashlib.sha256()
    for f in tests:
        th.update(posix(f).encode())''')

# 4. the rest of the old layout's names
edit('tools/cfd_kit.py', "from opencif import ENGINE  # noqa: E402,F401\n", "")
edit('tools/dissertation_figures.py', "    from opencif import renders as RR", "    from opencif.figures import proof as RR")
edit('provenance/rc7_traces/code/alignment_fixes.py', "    from opencif import renders as R", "    from opencif.figures import proof as R")
edit('tests/test_appendix_trace.py', "hashlib.md5(open(os.path.join(ROOT, 'opencif', 'reference', 'CIF_v2.py'), 'rb').read()).hexdigest()",
     "hashlib.md5(open(os.path.join(ROOT, 'reference', 'CIF_v2.py'), 'rb').read()).hexdigest()")
edit('tests/test_options_sheet.py', "open(os.path.join(ROOT, 'opencif', 'engine', 'layout.py')).read()",
     "open(os.path.join(ROOT, 'opencif', 'geometry', 'layout.py')).read()")
edit('tests/test_rc7.py', "m2src = open(os.path.join(ROOT, 'opencif', 'engine', 'measure2.py')).read()",
     "m2src = open(os.path.join(ROOT, 'opencif', 'checks', 'measure.py')).read()")
print('patch_v1: %d edits' % n_done)
