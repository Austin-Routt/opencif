"""opencif 1.0.0: the modules' own names in their docstrings and comments (the restructure moved and renamed the
files; the words naming them still said cifsched, measure2, renders, ...).  Only comments, docstrings and messages
change: every replacement below is checked to touch no code token.  Applied after patch_cli.py.
usage: python -I patch_names.py TREE"""
import io
import os
import re
import sys
import tokenize

T = sys.argv[1]
n_files, n_subs = 0, 0

FIRST = {   # module -> its new first docstring line prefix (old prefix -> new)
    'opencif/checks/measure.py': ('measure2.py -- ', 'opencif.checks.measure -- '),
    'opencif/published.py': ('devspecs.py -- ', 'opencif.published -- '),
    'opencif/physics/scenarios.py': ('robust.py -- ', 'opencif.physics.scenarios -- '),
    'opencif/physics/schedule.py': ('sched.py -- physics schedule for a CIF device with a rounded arms-meander head (clean-room proposal).',
                                    'opencif.physics.schedule -- the row schedule of a CIF device: the march at a constant lamina, the head and the exact closure.'),
    'opencif/figures/proof.py': ('renders.py -- ', 'opencif.figures.proof -- '),
    'opencif/figures/drawing.py': ('render.py -- ', 'opencif.figures.drawing -- '),
    'opencif/geometry/dxf.py': ('io_dxf.py -- ', 'opencif.geometry.dxf -- '),
    'opencif/geometry/shapes.py': ('geometry.py -- ', 'opencif.geometry.shapes -- '),
    'opencif/geometry/layout.py': ('layout.py -- ', 'opencif.geometry.layout -- '),
    'opencif/geometry/depicted.py': ('depicted.py -- ', 'opencif.geometry.depicted -- '),
    'opencif/geometry/dinh_inlet.py': ('dinh_inlet.py -- ', 'opencif.geometry.dinh_inlet -- '),
    'opencif/physics/network.py': ('network.py -- ', 'opencif.physics.network -- '),
    'opencif/api.py': ('api.py -- ', 'opencif.api -- '),
    'opencif/build.py': ('build.py -- ', 'opencif.build -- '),
    'opencif/pipeline.py': ('pipeline.py -- ', 'opencif.pipeline -- '),
    'opencif/spec.py': ('spec.py -- ', 'opencif.spec -- '),
    'opencif/compose.py': ('compose.py -- ', 'opencif.compose -- '),
}
SUBS = [   # (regex, replacement), applied to comments and string literals only
    (r'engine/cifsched\.py', 'opencif/physics/schedule.py'),
    (r'\bcifsched\b', 'schedule'),
    (r'\bmeasure2\b', 'measure'),
    (r'opencif/renders\.py', 'opencif/figures/proof.py'),
    (r'opencif\.renders\.', 'opencif.figures.proof.'),
    (r'\bgeometry\.corner_round_area\b', 'shapes.corner_round_area'),
    (r'\(geometry\.py\)', '(shapes.py)'),
    (r'\(layout\.py assembles', '(opencif.geometry.layout assembles'),
]
DEVSPECS_TEST = 'tests/test_appendix_trace.py'      # its alias for opencif.published, renamed with its messages


def sub_text(s):
    global n_subs
    for rx, rep in SUBS:
        s, k = re.subn(rx, rep, s)
        n_subs += k
    return s


for root, dirs, files in os.walk(T):
    dirs[:] = [d for d in dirs if d not in ('__pycache__', 'provenance', 'reference', 'foundation', 'reference_glossary')]
    for f in files:
        if not f.endswith('.py'):
            continue
        p = os.path.join(root, f)
        rel = os.path.relpath(p, T).replace(os.sep, '/')
        src0 = open(p, encoding='utf-8').read()
        src = src0
        if rel in FIRST:                                  # the first docstring line, before the word substitutions
            a, b = FIRST[rel]
            assert src.startswith('"""' + a), (rel, src[:80])
            src = '"""' + b + src[3 + len(a):]
            n_subs += 1
        # rewrite only COMMENT and STRING tokens, in place in the source text (from the last token back, so that the
        # offsets of the earlier ones stay valid); every other character is left as it is
        starts = [0]
        for ln in src.splitlines(keepends=True):
            starts.append(starts[-1] + len(ln))
        off = lambda rc: starts[rc[0] - 1] + rc[1]
        toks = [t for t in tokenize.generate_tokens(io.StringIO(src).readline) if t.type in (tokenize.COMMENT, tokenize.STRING)]
        new = src
        for t in reversed(toks):
            a_, b_ = off(t.start), off(t.end)
            assert new[a_:b_] == t.string, (rel, t)
            new = new[:a_] + sub_text(t.string) + new[b_:]
        if rel == DEVSPECS_TEST:
            k = new.count('devspecs')
            assert k == 3, k
            new = new.replace('devspecs', 'published')
            n_subs += k
        if new != src0:
            # the code tokens are the same (names, numbers, operators, in order)
            code = lambda s_: [(t.type, t.string) for t in tokenize.generate_tokens(io.StringIO(s_).readline)
                               if t.type not in (tokenize.COMMENT, tokenize.STRING, tokenize.NL, tokenize.NEWLINE,
                                                 tokenize.INDENT, tokenize.DEDENT)]
            if rel != DEVSPECS_TEST:
                assert code(src0) == code(new), rel
            open(p, 'w', encoding='utf-8').write(new)
            n_files += 1
print('patch_names: %d files, %d substitutions' % (n_files, n_subs))
