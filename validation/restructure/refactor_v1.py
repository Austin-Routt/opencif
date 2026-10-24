"""Restructure the frozen RC8 tree into the opencif 1.0.0 layout -- behaviour-preserving: files move, imports are
rewritten to package imports, the sys.path hacks go; no algorithm line changes.  Re-runnable: it always starts from
the RC8 extraction.   usage: python -I refactor_v1.py BASELINE_TREE OUT_TREE KERNELS_PY"""
import os
import re
import shutil
import sys

BASE, OUT, KERNELS = sys.argv[1:4]

MOVES = {   # old path (relative to the tree) -> new path
    'opencif/engine/cifsched.py': 'opencif/physics/schedule.py',
    'opencif/engine/network.py': 'opencif/physics/network.py',
    'opencif/engine/robust.py': 'opencif/physics/scenarios.py',
    'opencif/engine/geometry.py': 'opencif/geometry/shapes.py',
    'opencif/engine/layout.py': 'opencif/geometry/layout.py',
    'opencif/engine/dinh_inlet.py': 'opencif/geometry/dinh_inlet.py',
    'opencif/engine/depicted.py': 'opencif/geometry/depicted.py',
    'opencif/engine/io_dxf.py': 'opencif/geometry/dxf.py',
    'opencif/engine/devspecs.py': 'opencif/published.py',
    'opencif/engine/render.py': 'opencif/figures/drawing.py',
    'opencif/renders.py': 'opencif/figures/proof.py',
    'opencif/engine/measure2.py': 'opencif/checks/measure.py',
    'opencif/engine/drc_audit.py': 'opencif/checks/drc.py',
    'opencif/reference/CIF_v2.py': 'reference/CIF_v2.py',
}
MODMAP = {   # old top-level module name -> (package, module)
    'cifsched': ('opencif.physics', 'schedule'),
    'network': ('opencif.physics', 'network'),
    'robust': ('opencif.physics', 'scenarios'),
    'geometry': ('opencif.geometry', 'shapes'),
    'layout': ('opencif.geometry', 'layout'),
    'dinh_inlet': ('opencif.geometry', 'dinh_inlet'),
    'depicted': ('opencif.geometry', 'depicted'),
    'io_dxf': ('opencif.geometry', 'dxf'),
    'devspecs': ('opencif', 'published'),
    'render': ('opencif.figures', 'drawing'),
    'measure2': ('opencif.checks', 'measure'),
}
PKG_INIT = {
    'opencif/physics/__init__.py': '"""The physics: the kernels (Appendix A), the row schedule (march and closure), the 1-D flow network, and the\n'
                                   'model scenarios."""\n',
    'opencif/geometry/__init__.py': '"""Drawing a device: the shapes of its parts, the layout of legs, heads, turns, inlets and outlets, the\n'
                                    'published devices\' depicted parts, and the DXF writer."""\n',
    'opencif/checks/__init__.py': '"""Checking a drawn device: the measurement of the drawn polygon, the design-rule check of the DXF, the\n'
                                  'physics auditor, and the design code\'s reports and acceptance validators."""\n',
    'opencif/figures/__init__.py': '"""Figures: each build\'s proof figures and the drawing styles they use."""\n',
}

if os.path.exists(OUT):
    shutil.rmtree(OUT)
shutil.copytree(BASE, OUT, ignore=shutil.ignore_patterns('__pycache__'))
for old, new in MOVES.items():
    os.makedirs(os.path.dirname(os.path.join(OUT, new)), exist_ok=True)
    shutil.move(os.path.join(OUT, old), os.path.join(OUT, new))
for d in ('opencif/engine', 'opencif/reference'):
    left = os.listdir(os.path.join(OUT, d))
    assert not [f for f in left if f.endswith('.py')], (d, left)
    shutil.rmtree(os.path.join(OUT, d))
for p, txt in PKG_INIT.items():
    with open(os.path.join(OUT, p), 'w') as fh:
        fh.write(txt)
shutil.copy(KERNELS, os.path.join(OUT, 'opencif/physics/kernels.py'))

IMPORT_RE = re.compile(r'^(?P<ind>\s*)import (?P<body>[A-Za-z0-9_]+(?: as [A-Za-z0-9_]+)?(?:\s*,\s*[A-Za-z0-9_]+(?: as [A-Za-z0-9_]+)?)*)'
                       r'(?P<rest>\s*(?:#.*)?)$')
SPECIAL = [   # exact line rewrites (stripped line -> replacement text without indentation); None deletes the line
    ("from . import ENGINE  # noqa: F401  (puts the engine on sys.path)", None),
    ("from . import ENGINE  # noqa: F401", None),
    ("from . import ENGINE  # noqa", None),
    ("from . import renders as RN", "from .figures import proof as RN"),
    ("from opencif.engine import cifsched as S        # devicekit.use_rc6() puts the package's rc6 first on the path",
     "from opencif.physics import schedule as S"),
    ("sys.path.insert(0, os.path.join(PKG, 'reference'))                      # CIF_v2.py v2.1.0, byte-frozen: kernels only",
     None),
    ("import CIF_v2 as C          # noqa: E402  official kernels",
     "from opencif.physics import kernels as C    # the kernels, bit-identical to reference/CIF_v2.py (tests/test_kernels.py)"),
    ("from opencif.reference import CIF_v2 as C",
     "sys.path.insert(0, os.path.join(ROOT, 'reference'))\nimport CIF_v2 as C                      # noqa: E402  the frozen design code"),
]
SPECIAL_MAP = {a: b for a, b in SPECIAL}
used_special = set()
changes = {}


def rewrite_import(m):
    ind, body, rest = m.group('ind'), m.group('body'), m.group('rest') or ''
    parts = [p.strip() for p in body.split(',')]
    names = []
    for p in parts:
        nm, _, al = p.partition(' as ')
        names.append((nm.strip(), al.strip() or None))
    if not any(nm in MODMAP for nm, _ in names):
        return None
    keep = [(nm, al) for nm, al in names if nm not in MODMAP]
    out = []
    if keep:
        out.append(ind + 'import ' + ', '.join(nm + (' as ' + al if al else '') for nm, al in keep))
    for nm, al in names:
        if nm in MODMAP:
            pkg, mod = MODMAP[nm]
            out.append('%sfrom %s import %s as %s' % (ind, pkg, mod, al or nm))
    if rest.strip():
        out[0] = out[0] + '  ' + rest.strip()
    return out


ENGINE_PATH_RE = [   # sys.path entries that pointed into the package or its engine
    re.compile(r",\s*os\.path\.join\((\w+), 'opencif', 'engine'\)"),
    re.compile(r",\s*os\.path\.join\((\w+), 'opencif'\)(?=\s*[,\]\)])"),
]
for root, dirs, files in os.walk(OUT):
    dirs[:] = [d for d in dirs if d != '__pycache__']
    for f in files:
        if not f.endswith('.py'):
            continue
        p = os.path.join(root, f)
        rel = os.path.relpath(p, OUT)
        if rel == 'reference/CIF_v2.py' or rel == 'opencif/physics/kernels.py':
            continue
        src = open(p, encoding='utf-8').read()
        lines = src.split('\n')
        new = []
        n_ch = 0
        for ln in lines:
            s = ln.strip()
            if s in SPECIAL_MAP:
                used_special.add(s)
                rep = SPECIAL_MAP[s]
                n_ch += 1
                if rep is None:
                    continue
                ind = ln[:len(ln) - len(ln.lstrip())]
                new.extend(ind + x for x in rep.split('\n'))
                continue
            m = IMPORT_RE.match(ln)
            if m:
                r = rewrite_import(m)
                if r is not None:
                    new.extend(r)
                    n_ch += 1
                    continue
            ln2 = ln
            if 'sys.path' in ln or 'for p in (' in ln:
                for rx in ENGINE_PATH_RE:
                    ln2 = rx.sub('', ln2)
            if ln2 != ln:
                n_ch += 1
            new.append(ln2)
        out = '\n'.join(new)
        if out != src:
            open(p, 'w', encoding='utf-8').write(out)
            changes[rel] = n_ch
missing = [a for a, _ in SPECIAL if a not in used_special]
print('rewrote %d files (%d edits); special rules unused: %s' % (len(changes), sum(changes.values()), missing))
for k in sorted(changes):
    print('  %3d  %s' % (changes[k], k))
