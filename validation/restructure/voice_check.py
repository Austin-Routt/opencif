"""voice_check.py: prove that a voiced copy of a module changes only its comments and docstrings, and scan its prose.

Usage: python -I voice_check.py BASE_DIR NEW_DIR [relpath ...]

For every .py file under NEW_DIR (or the relpaths given), the check parses the file in BASE_DIR and in NEW_DIR,
removes every docstring, and requires the two abstract syntax trees to be identical.  It also
requires that no string constant other than a docstring differs, that the test-anchored source passages are still
present, and that the file compiles.  The prose scan lists the constructions that the author's style sheet excludes,
so that a reviewer can read each one.  The exit code is 0 when every file passes the AST check and 1 otherwise.
"""
import ast, io, re, sys, tokenize, pathlib

# Source passages that the acceptance suites read from the package source; a voiced file must keep each one.
ANCHORS = {
    'opencif/build.py': ['kernel_pair=KERNEL_PAIR'],
    'opencif/pipeline.py': ['kernel_pair=BU.KERNEL_PAIR', 'status_codes=codes', "status['source_trace']",
                            'paper.targets_fitted: %s'],
    'opencif/api.py': ['kernel_pair=BU.KERNEL_PAIR'],
    'opencif/checks/measure.py': ['Lm = (area - a_round) / wm + excl'],
    'opencif/geometry/layout.py': ['clear=150.0'],
    'opencif/__main__.py': ['_print_status(rec)'],
    'opencif/checks/auditor.py': ['spec.F_target - 1e-3', 'abs(F_em - F_march) < 0.02', 'max(8e-3, 4.0 * fbar)',
                                  'mismatch=3.0', 'Le = 0.05 * Re * Dh', 'Le <= 0.5 * geom["Lambda"]',
                                  'ubar = Q_SI / (wc0_SI * h_SI)'],
}

BRITISH = r'\b(centre[sd]?|centring|realis\w*|modell\w*|licence\w*|neighbour\w*|grey\w*|metres?|behaviour\w*|' \
          r'colour\w*|favour\w*|normalis\w*|minimis\w*|maximis\w*|organis\w*|recognis\w*|parameteris\w*|analys(e|ed|ing)|' \
          r'catalogue\w*|programme\w*|labelled|labelling|travell\w*|fibre\w*|optimis\w*|summaris\w*|generalis\w*|' \
          r'characteris\w*|equalis\w*|stabilis\w*|utilis\w*|emphasis(e|ed|es)|prioritis\w*|visualis\w*|initialis\w*|' \
          r'finalis\w*|serialis\w*|digitis\w*|discretis\w*|regularis\w*|symmetris\w*|tonne\w*|aluminium|judgement|' \
          r'artefact\w*|whilst|amongst)\b'
EXCLUDED = [
    (r'\u2014|–', 'em or en dash'),
    (r'(?<!-)\s--\s(?!-)', 'double-hyphen dash'),
    (r'\byou\b|\byour\b|\byours\b', '"you"'),
    (r'\b(notably|importantly|crucially|essentially|straightforward\w*|simply|clearly|robust\w*|ensures?|ensuring|'
     r'leverag\w*|comprehensive|overall|interestingly|specifically|namely|in order to|it is worth noting)\b',
     'excluded word'),
    (r"\b\w+n't\b|\b(it's|that's|there's|we're|we've|let's)\b", 'contraction'),
    (r'\*\*', 'bold'),
    (r'\bnot\b[^.;]{1,60}\bbut\b', '"not X but Y"'),
]


def _strip_docstrings(tree):
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = node.body
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                del body[0]                                  # a docstring added or removed is prose too
    return tree


def _docstrings(tree):
    out = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            d = ast.get_docstring(node, clean=False)
            if d:
                out.append((node.body[0].lineno, d))
    return out


def _code_tokens(src):
    """Return the code's tokens with their columns, leaving out comments, blank lines, and docstrings.

    The AST comparison does not see a change of spacing inside a line of code, so this second comparison holds every
    code token to its column as well."""
    tree = ast.parse(src)
    doc_lines = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) \
                    and isinstance(first.value.value, str):
                doc_lines.add(first.lineno)
    out = []
    skip_newline = False
    for t in tokenize.generate_tokens(io.StringIO(src).readline):
        if t.type in (tokenize.COMMENT, tokenize.NL, tokenize.INDENT, tokenize.DEDENT, tokenize.ENCODING,
                      tokenize.ENDMARKER):
            continue
        if t.type == tokenize.STRING and t.start[0] in doc_lines:
            skip_newline = True
            continue
        if t.type == tokenize.NEWLINE and skip_newline:
            skip_newline = False
            continue
        skip_newline = False
        # a NEWLINE token's column follows any trailing comment, so only its place in the sequence is compared
        out.append((t.type, t.string, -1 if t.type == tokenize.NEWLINE else t.start[1]))
    return out


def _comments(src):
    return [(t.start[0], t.string) for t in tokenize.generate_tokens(io.StringIO(src).readline)
            if t.type == tokenize.COMMENT]


def _prose_hits(text):
    """Return the excluded constructions in one comment or docstring, with code spans and identifiers removed."""
    t = re.sub(r'`[^`]*`', ' ', text)                       # code spans
    t = re.sub(r'\b\w*_\w*\b', ' ', t)                       # identifiers with underscores
    t = re.sub(r'\b[\w.]+\(', ' ', t)                        # calls
    t = re.sub(r"'[^'\n]*'", ' ', t)                         # quoted keys and messages
    hits = [why for pat, why in EXCLUDED if re.search(pat, t, re.I if why != 'bold' else 0)]
    if re.search(BRITISH, t, re.I):
        hits.append('British spelling: ' + re.search(BRITISH, t, re.I).group(0))
    return hits


def check(base, new, rel):
    a_src = (base / rel).read_text(encoding='utf-8')
    b_src = (new / rel).read_text(encoding='utf-8')
    problems = []
    try:
        compile(b_src, rel, 'exec')
    except SyntaxError as e:
        return ['does not compile: %s' % e], []
    a_tree, b_tree = ast.parse(a_src), ast.parse(b_src)
    if ast.dump(_strip_docstrings(a_tree)) != ast.dump(_strip_docstrings(b_tree)):
        problems.append('AST differs outside the docstrings')
    elif _code_tokens(a_src) != _code_tokens(b_src):
        ta, tb = _code_tokens(a_src), _code_tokens(b_src)
        i = next((k for k, (x, y) in enumerate(zip(ta, tb)) if x != y), min(len(ta), len(tb)))
        problems.append('the layout of the code changed near %r' % (tb[i][1] if i < len(tb) else 'the end'))
    for anchor in ANCHORS.get(rel, []):
        if anchor in a_src and anchor not in b_src:
            problems.append('test anchor lost: %r' % anchor)
    if rel == 'opencif/pipeline.py':
        i = b_src.index('def run(user, outdir=None, verbose=True, render=True):')
        if 'S.clear_caches()' not in b_src[i:i + 300]:
            problems.append('test_rc7 BUILD-1: S.clear_caches() is no longer within 300 characters of def run(')
    prose = []
    for ln, c in _comments(b_src):
        for h in _prose_hits(c):
            prose.append((ln, h, c.strip()[:110]))
    for ln, d in _docstrings(ast.parse(b_src)):
        for k, line in enumerate(d.splitlines()):
            for h in _prose_hits(line):
                prose.append((ln + k, h, line.strip()[:110]))
    return problems, prose


def main(argv):
    base, new = pathlib.Path(argv[1]), pathlib.Path(argv[2])
    rels = argv[3:] or sorted(str(p.relative_to(new)) for p in new.rglob('*.py'))
    bad = 0
    for rel in rels:
        problems, prose = check(base, new, rel)
        print('%s: %s; %d prose flags' % (rel, 'AST equal' if not problems else '; '.join(problems), len(prose)))
        for ln, why, text in prose:
            print('    %5d  %-28s %s' % (ln, why, text))
        bad += bool(problems)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
