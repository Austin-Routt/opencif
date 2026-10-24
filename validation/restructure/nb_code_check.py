"""nb_code_check.py: compare the code cells of two versions of the tutorial notebook.

Usage: python -I nb_code_check.py OLD.ipynb NEW.ipynb

Each code cell is parsed, its comments are dropped by the parser, and every string constant is replaced by a
placeholder, so two cells match when they differ only in comments, docstrings, and the text of their strings. The
cells are aligned in order, and the report lists every new cell, every removed cell, every cell whose code changed,
and every changed string whose count of %-conversions differs (a changed format string would fail at run time).
"""
import ast, difflib, json, re, sys


class _Blank(ast.NodeTransformer):
    def visit_Constant(self, node):
        if isinstance(node.value, str):
            return ast.copy_location(ast.Constant('S'), node)
        return node

    def visit_JoinedStr(self, node):
        return ast.copy_location(ast.Constant('S'), node)


def _cells(path):
    nb = json.load(open(path, encoding='utf-8'))
    out = []
    for c in nb['cells']:
        if c['cell_type'] != 'code':
            continue
        src = ''.join(c['source'])
        tree = ast.parse(src)
        strings = [n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)]
        out.append((ast.dump(_Blank().visit(tree)), strings, src))
    return out


def _conv(s):
    return len(re.findall(r'%(?!%)[-+ #0]*\d*(?:\.\d+)?[sdifgreEx]', s))


def main(argv):
    a, b = _cells(argv[1]), _cells(argv[2])
    sm = difflib.SequenceMatcher(a=[x[0] for x in a], b=[x[0] for x in b], autojunk=False)
    bad = 0
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op == 'equal':
            for i, j in zip(range(i1, i2), range(j1, j2)):
                for s_old, s_new in zip(a[i][1], b[j][1]):
                    if s_old != s_new and _conv(s_old) != _conv(s_new):
                        print('cell %d -> %d: a format string lost or gained a conversion: %r -> %r'
                              % (i, j, s_old[:60], s_new[:60]))
                        bad += 1
            continue
        for i in range(i1, i2):
            print('old code cell %d changed or removed: %s' % (i, a[i][2].strip().splitlines()[0][:100]))
        for j in range(j1, j2):
            print('new code cell %d added or changed: %s' % (j, b[j][2].strip().splitlines()[0][:100]))
        bad += op != 'insert'
    print('%d old and %d new code cells; %d problems' % (len(a), len(b), bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
