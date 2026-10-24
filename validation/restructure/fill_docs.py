"""fill_docs.py: fill the placeholders of the release's documents once the evidence exists.

Usage: python -I fill_docs.py TREE VALUES_JSON

The documents carry placeholders such as {{STAMP}} and {{AUDIT_SUMMARY}} wherever a number or a sentence depends on
the final evidence. The documents are README.md, CITATION.cff, the workflow that checks the code stamp
(.github/workflows/tests.yml), every docs/*.md, and the errata list of the frozen texts
(foundation/ERRATA_NEXT_UNFREEZE.md, whose proposed text names the code stamp). The code stamp, the tests stamp,
the Python version, and the line counts of the module map ({{LINES:path,...}}) are computed here from the tree; every
other placeholder takes its text from VALUES_JSON, which make_values_v1.py writes from validation/RELEASE.json, the
suites' logs, and the audits. A placeholder without a value stops the script, so no document is published with one
left in it.
"""
import json
import os
import re
import sys

TREE, VALUES = sys.argv[1], sys.argv[2]
sys.path.insert(0, TREE)
from opencif import pipeline as PL                                   # noqa: E402

values = json.load(open(VALUES, encoding='utf-8'))
env = PL.environment_stamp() if hasattr(PL, 'environment_stamp') else PL.ENV
values.setdefault('STAMP', PL.CODE['sha256'])
values.setdefault('TESTS_STAMP', env['tests_sha256'])
values.setdefault('PYTHON', env['python'])


def n_lines(rel):
    with open(os.path.join(TREE, 'opencif', *rel.split('/')), encoding='utf-8') as fh:
        return sum(1 for _ in fh)


def lines_value(spec):
    return ', '.join('{:,}'.format(n_lines(r.strip())) for r in spec.split(','))


targets = ['README.md', 'CITATION.cff', '.github/workflows/tests.yml', 'foundation/ERRATA_NEXT_UNFREEZE.md'] + sorted('docs/' + f for f in os.listdir(os.path.join(TREE, 'docs'))
                                                 if f.endswith('.md'))
missing, n_filled = set(), 0
for rel in targets:
    path = os.path.join(TREE, *rel.split('/'))
    if not os.path.exists(path):
        continue
    text = open(path, encoding='utf-8').read()

    def sub(m):
        global n_filled
        key = m.group(1)
        if key.startswith('LINES:'):
            n_filled += 1
            return lines_value(key[len('LINES:'):])
        if key in values:
            n_filled += 1
            return str(values[key])
        missing.add(key)
        return m.group(0)
    new = re.sub(r'\{\{([A-Za-z0-9_:/.,]+)\}\}', sub, text)
    if new != text:
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(new)
if missing:
    raise SystemExit('fill_docs: no value for %s' % ', '.join(sorted(missing)))
print('fill_docs: %d placeholders filled in %d documents' % (n_filled, len(targets)))
