"""patch_voice.py: put the comments and docstrings of opencif 1.0.0 in the author's voice.

Usage: python -I patch_voice.py TREE VOICE_DIR [SUBDIR ...]   (default SUBDIR: opencif)

VOICE_DIR holds a voiced copy of each module (VOICE_DIR/opencif/<path>.py) and of each test suite
(VOICE_DIR/tests/<name>.py). Each copy replaces the tree's module only
when the two modules have the same abstract syntax tree once their docstrings are removed, so the rewrite can change
comments and docstrings and nothing else (voice_check.py states the rule). The modules whose sources live in tools_v1
(the kernels, the auditor and its ruler, the design code's reports, and the teaching module) are voiced at their
source and are not in VOICE_DIR. The package is voiced after patch_lint.py, and the tests after patch_tests.py.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import voice_check as VC                                          # noqa: E402

tree, vdir = sys.argv[1], sys.argv[2]
subdirs = sys.argv[3:] or ['opencif']
done, bad = 0, []
for root, _, files in (w for sub in subdirs for w in os.walk(os.path.join(vdir, sub))):
    for fn in sorted(files):
        if not fn.endswith('.py'):
            continue
        rel = os.path.relpath(os.path.join(root, fn), vdir).replace(os.sep, '/')
        problems, _ = VC.check(VC.pathlib.Path(tree), VC.pathlib.Path(vdir), rel)
        if problems:
            bad.append((rel, problems))
            continue
        with open(os.path.join(vdir, rel), 'rb') as fh:
            data = fh.read()
        with open(os.path.join(tree, rel), 'wb') as fh:
            fh.write(data)
        done += 1
if bad:
    raise SystemExit('patch_voice: %d module(s) would change code: %s' % (len(bad), bad))
print('patch_voice: %d files of %s voiced, each proven to change comments and docstrings only' % (done, ', '.join(subdirs)))
