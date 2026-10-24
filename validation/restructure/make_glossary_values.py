"""Write reference/glossary_values.json: the numeric macros of the v114 glossary's preamble that tests/test_appendix_trace.py
pins the code to, so that a copy of the repository without the foundation texts runs the same checks.  With the
preamble present (reference_glossary/PREAMBLE_FINAL.tex, as the repository ships it), the test requires its md5 and
that this file holds exactly its values.   usage: python -I make_glossary_values.py TREE"""
import hashlib
import json
import os
import re
import sys

TREE = sys.argv[1]
USED = ['FclosedVone', 'FclosedVthree', 'dcritGGIn', 'dinhBandHi', 'dinhBandLo', 'dinhCollapseLo', 'dinhKappaLocal',
        'footUwc', 'hillKmax', 'hillN', 'hillUhalf', 'kappCIhalfPct', 'kappCIhi', 'kappCIlo', 'kappOverRootThree',
        'kappVal', 'plateauUwc', 'squareDuctR']
pre_p = os.path.join(TREE, 'reference_glossary', 'PREAMBLE_FINAL.tex')
raw = open(pre_p, 'rb').read()
pre = raw.decode('utf-8')
mac = {m.group(1): m.group(2) for m in re.finditer(r'\\newcommand\{\\(\w+)\}\{([0-9.]+)\}', pre)}
cap = re.search(r'\\newcommand\{\\fgapDesignCap\}\{([0-9.]+)\\times10\^\{(-?\d+)\}\}', pre)
missing = [k for k in USED if k not in mac]
if missing or cap is None:
    raise SystemExit('make_glossary_values: not in the preamble: %s' % (missing or ['fgapDesignCap']))
out = dict(source='the numeric macros of the dissertation glossary\'s preamble (v114, PREAMBLE_FINAL.tex, md5 %s) that '
                  'tests/test_appendix_trace.py pins the code to; the preamble itself is foundation/PREAMBLE_FINAL.tex'
                  % hashlib.md5(raw).hexdigest(),
           preamble_md5=hashlib.md5(raw).hexdigest(),
           macros={k: mac[k] for k in USED},
           fgapDesignCap=[cap.group(1), cap.group(2)])
p = os.path.join(TREE, 'reference', 'glossary_values.json')
with open(p, 'w', encoding='utf-8') as fh:
    json.dump(out, fh, indent=1, sort_keys=True)
    fh.write('\n')
print('wrote %s: %d macros' % (p, len(USED) + 1))
