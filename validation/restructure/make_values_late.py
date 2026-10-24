"""make_values_late.py: write the three values of the release's documents that wait for work done after the build
(docs/VALIDATION.md sections 6 and 7): the device packages' results, the executed tutorial's, and the independent
reviews'. fill_docs.py fills them in a second pass, after the values of make_values_v1.py.

Usage: python -I make_values_late.py TRACES_OUT TUTORIAL_IPYNB OUT.json [SPOT_CHECK_LOG]

TRACES_OUT      the folder into which tools/device_package.py ran the four device packages (opencif_<device>/out/
                run_log.txt in each), as the traces notebook does
TUTORIAL_IPYNB  the executed tutorial (notebooks/opencif_tutorial.ipynb after its full run)
SPOT_CHECK_LOG  the log of the spot rebuild and the new design of docs/EVALUATION.md sections 4 and 5, run as printed
                (criterion 6 of its section 8), with one "RESULT <name>: ..." line per rebuilt record and one for
                section 5

Every number of the first two values is read from these files, and the run stops if a package did not finish, if a
photomask is not byte-identical to the release's, or if a cell of the tutorial is in error. The third value describes
the reviews (their reports and the record of what was done are in the release's evidence assets).
"""
import datetime
import json
import os
import re
import sys

TRACES, NB, OUT = sys.argv[1:4]
SPOT = sys.argv[4] if len(sys.argv) > 4 else None
DEVICES = ['gifford2014_fig6', 'gifford2014_fig4', 'dinh2024_v1', 'xia2016_rs7']
TS = re.compile(r'^(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)  (.*)$')


def stamp(s):
    return datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S')


def device(dev):
    log = os.path.join(TRACES, 'opencif_' + dev, 'out', 'run_log.txt')
    lines = [TS.match(ln.rstrip('\n')) for ln in open(log, encoding='utf-8')]
    lines = [(stamp(m.group(1)), m.group(2)) for m in lines if m]
    done = [t for _, t in lines if t.startswith('done: ')]
    if not done:
        raise SystemExit('%s: the package did not finish (no "done:" line in %s)' % (dev, log))
    steps = done[-1][len('done: '):]
    masks = [t for _, t in lines if t.startswith('build: photomask ')]
    same = [t for t in masks if 'byte-identical to the release' in t]
    if not masks or len(same) != len(masks):
        raise SystemExit('%s: %d of %d photomasks byte-identical' % (dev, len(same), len(masks)))
    regs = [t for _, t in lines if t.startswith('trace: every value of RC8')]
    stamp_ok = any('code stamp 9e4eb55a0affeef3' in t and 'the release' in t for _, t in lines)
    vers_ok = any(t.startswith('every package at the version') for _, t in lines)
    if not (stamp_ok and vers_ok):
        raise SystemExit('%s: the check step did not confirm the code stamp and the tested versions' % dev)
    extra = []
    for _, t in lines:
        m = re.search(r'(\d+) checks within their regression limits', t)
        if m:
            extra.append('%s measurements within their regression limits' % m.group(1))
        if t.startswith('taper: every V13 result equal'):
            extra.append('every result of the taper-direction test (V13) equal to RC8\'s, with the nested reading\'s '
                         'neighbouring leg PASS, reading B\'s FAIL and every leg-4 control PASS')
        if 'claim (leg 4) lies on the inset within its regression limits' in t:
            extra.append('the claim\'s leg on the inset within its regression limits')
    minutes = (lines[-1][0] - lines[0][0]).total_seconds() / 60
    return {'dev': dev, 'steps': steps, 'masks': len(masks), 'regs': len(regs), 'extra': extra,
            'minutes': minutes, 'date': lines[0][0].strftime('%Y-%m-%d')}


def words(n):
    return {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six'}.get(n, str(n))


D = [device(d) for d in DEVICES]
parts = []
for d in D:
    s = '`%s` (%.0f min): %s photomask%s byte-identical' % (d['dev'], d['minutes'], words(d['masks']),
                                                         '' if d['masks'] == 1 else 's')
    if d['regs']:
        s += ', %s registration record%s equal to RC8\'s in every value' % (words(d['regs']),
                                                                          '' if d['regs'] == 1 else 's')
    for e in d['extra']:
        s += ', ' + e
    parts.append(s)
n_masks = sum(d['masks'] for d in D)
V = {}
V['DEVICES_RESULT'] = (
    'The four device packages were assembled from this tree by `tools/device_package.py` and run with the '
    'publishers\' figures on %s, each in its own process. Each first checked the tested versions, the code stamp and '
    'the sha256 of every input, then built its photomasks, registered them on the papers\' micrographs and rendered '
    'the overlays, and every one ended its last step. All %d photomasks are byte-identical to the release\'s, which '
    'are RC8\'s. Per package: %s. The packages\' outputs, which embed the publishers\' figures, are among the '
    'release\'s evidence assets rather than in this repository. The first run of the Xia package stopped at its trace '
    'step on an error of the restructure in five trace scripts, which was corrected (CHANGELOG item 14) before all '
    'four packages were run again; the results above are from that second run. After it, the kits\' source notes '
    '(`inputs/SOURCE.md`) were updated to say that the figures are distributed in `notebooks/figures/`, which changed '
    'the sha256 of each kit\'s source note in its `MANIFEST.json` and `reference/EXPECT.json` and no other entry, and '
    'each kit then assembled from the published tree with the figures of `notebooks/figures/` and passed its check step.'
    % (D[0]['date'], n_masks, '; '.join(parts)))

# ---- the executed tutorial
nb = json.load(open(NB, encoding='utf-8'))
code = [c for c in nb['cells'] if c['cell_type'] == 'code']
errors = [c for c in code for o in c.get('outputs', []) if o.get('output_type') == 'error']
unrun = [c for c in code if c.get('execution_count') is None]
if errors or unrun:
    raise SystemExit('tutorial: %d cells in error and %d not run' % (len(errors), len(unrun)))
t = []
for c in code:
    ex = c.get('metadata', {}).get('execution', {})
    for k in ('iopub.execute_input', 'shell.execute_reply'):
        if k in ex:
            t.append(datetime.datetime.fromisoformat(ex[k].replace('Z', '+00:00')))
minutes = (max(t) - min(t)).total_seconds() / 60 if t else float('nan')


def text_of(cell):
    out = []
    for o in cell.get('outputs', []):
        if o.get('output_type') == 'stream':
            out.append(''.join(o['text']))
        elif 'data' in o and 'text/plain' in o['data']:
            out.append(''.join(o['data']['text/plain']))
    return ''.join(out)


def cell_with(src_part):
    hits = [c for c in code if src_part in ''.join(c['source'])]
    if len(hits) != 1:
        raise SystemExit('tutorial: %d cells contain %r' % (len(hits), src_part))
    return text_of(hits[0])


repro = cell_with("'equal to the release\\'s'")
rows = [ln for ln in repro.splitlines()[1:] if ln.strip() and not ln.startswith(' ')]
equal = [ln for ln in rows if ln.rstrip().endswith('yes, byte for byte')]
if not rows or len(equal) != len(rows):
    raise SystemExit('tutorial: %d of %d rebuilt photomasks equal to the release\'s:\n%s' % (len(equal), len(rows),
                                                                                            repro))
scratch = cell_with('again = AU.audit_run(first.out')
m = re.search(r"standalone audit (\S+) in \d+ s; every check's verdict equal to that of the build's own audit: (\w+)",
              scratch)
if not m or m.group(2) != 'True':
    raise SystemExit('tutorial: the audit from scratch did not give the same verdicts:\n' + scratch)
differs = re.search(r'checks whose recorded value differs [^:]*: (.*)', scratch).group(1).strip()
first = cell_with("first = opencif.generate(")
g = (re.search(r'verdict  : (\S+)', first), re.search(r'exit code: (\d+)', first),
     re.search(r'physics audit of the photomask file: (\S+)', first))
if not all(g) or g[1].group(1) != '0':
    raise SystemExit('tutorial: the new design of section 6 did not pass:\n' + first)
gen = 'verdict %s, exit %s, physics audit %s' % (g[0].group(1), g[1].group(1), g[2].group(1))
raw = open(NB, encoding='utf-8').read()
paths = sorted(set(re.findall(r'(/home/|/tmp/|/Users/|/root/|[A-Z]:\\\\)', raw)))
if paths:
    raise SystemExit('tutorial: absolute paths in the notebook: %s' % paths)
V['TUTORIAL_RESULT'] = (
    'The tutorial (`notebooks/opencif_tutorial.ipynb`) was executed in full (`QUICK = False`) on this tree with the '
    'tested versions in %.0f min, and none of its %d code cells is in error. The published notebook is that executed '
    'one, with the text of its section 12 updated afterwards to say that the figures are distributed (no code cell '
    'or output changed), and it contains no absolute path of the machine that ran it. Its section 6 built a new device from D = 4 µm '
    'and F = 0.5 (%s). Its section 7.1 audited that photomask again from scratch, and the audit gave %s, with every '
    'check\'s verdict equal to that of the build\'s own audit (the checks whose recorded value differs in its last '
    'digits or in the row where it peaks: %s). Its section 10 rebuilt the %s published devices, and every photomask '
    'is byte-identical to the release\'s.' % (minutes, len(code), gen, m.group(1), differs, words(len(rows))))

V['REVIEW_RESULT'] = (
    'Two reviewers worked independently on the tree as published and on the build\'s evidence (2026-10-23): one on the '
    'physics, the drawing and its measurement, and the gates, and the other on reproducibility, the documents against '
    'the records, the tutorial, and the repository\'s hygiene. Both found the physics sound. Every finding was checked '
    'against the code or the records before it was acted on. The release keeps its code stamp, because its '
    'photomasks, records and audits were built and validated at that stamp, so a finding about the package is answered '
    'by documents that say what the code does, and the code change is listed for the next release '
    '(`docs/LIMITATIONS.md` item 31). The documents, the tests\' workflow, the tools outside the package and the '
    'provenance were corrected where needed (CHANGELOG item 14). A third reviewer re-audited every row of '
    '`docs/ALIGNMENT.md` against `foundation/` and the restructured code (criterion 5 of `docs/EVALUATION.md` section '
    '8). It found one row to relabel, AA16, whose floor start departed from spec v5 §10 and M6 (row S3, erratum D108), '
    'and ten smaller corrections; after they were made, its re-check found every row correctly labelled. At the '
    'author\'s direction the specification was then corrected in place (revision 5.1, `foundation/README.md`). The '
    'two suites that read the documents were run again after every change, with identical logs. The reviews and the '
    'record of what was done with each finding are among the release\'s evidence assets.')

if SPOT:
    res = [ln.strip() for ln in open(SPOT, encoding='utf-8') if ln.startswith('RESULT ')]
    rebuilt = [r for r in res if not r.startswith('RESULT section 5')]
    sec5 = [r for r in res if r.startswith('RESULT section 5')]
    if (len(rebuilt) != 2 or not all(r.endswith('every compared field equal') for r in rebuilt) or len(sec5) != 1
            or not sec5[0].endswith('as the guide says')):
        raise SystemExit('the spot check did not pass:\n' + '\n'.join(res))
    names = ' and '.join('`%s`' % r.split()[1].rstrip(':') for r in rebuilt)
    V['REVIEW_RESULT'] += (
        ' With these, the release meets the six criteria of `docs/EVALUATION.md` section 8: sections 1 to 3 of this '
        'page and `docs/RESULTS.md` show criteria 1 to 4, the re-audit shows criterion 5, and for criterion 6 the commands of its sections 4 '
        'and 5 were run as printed: the rebuilt %s equal their records in every compared field (the verdict, the '
        'realised cutoff, the network, the physics audit, every check\'s value, and the photomask\'s sha256), and the '
        'new design passes and the three refusals refuse for the reasons printed (the comments on two of them were '
        'corrected to name the refusal on w_min at row 1 that comes first).' % names)

json.dump(V, open(OUT, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
for k, v in V.items():
    print('%s: %s\n' % (k, v))
