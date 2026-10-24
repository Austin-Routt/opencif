"""make_validation.py: write the release's evidence summary, validation/ of the opencif 1.0.0 tree.

Usage: python -I make_validation.py TREE RC8_FROZEN_JSON PLANS_RC8_JSON PLANS_V1_JSON RECORDS_JSON SUITES_DIR
                                    AUDITS_DIR TOOLS_DIR

TREE             the final tree, whose runs/ holds the rebuild (runs/rebuild_summary.jsonl and one folder per build)
RC8_FROZEN_JSON  the RC8 photomasks (name -> sha256, exit code, N, F), rc8ref/frozen_photomasks.json
PLANS_*_JSON     the plan fingerprints of every shipped specification in RC8's tree and in this one (regress_plans.py)
RECORDS_JSON     the comparison of every record with RC8's (compare_records.py)
SUITES_DIR       the logs of the acceptance suites run on this tree (<suite>.log, and summary.txt)
AUDITS_DIR       the standalone audits (<name>.json, from audit_run on each finished build)
TOOLS_DIR        the scripts that made this tree from RC8's frozen tree (tools_v1)
RC8_RUNS         RC8's runs/ (the two arrays' RC8 records, whose photomask sha256 the frozen list does not carry)

It writes validation/RELEASE.json (the stamps, the counts, and one entry per photomask), validation/audit_table.csv
(every check of every audit, in the build and from scratch), validation/photomasks_sha256.txt,
validation/rc8_comparison/ (its four comparison inputs, copied), validation/logs/ (the suites' logs, with the tree's
absolute folder written <repo>, and the rebuild's summary), and validation/restructure/ (the scripts). Nothing is computed here that a log or a record does not already hold, so every number in RELEASE.json
can be traced to its source file.
"""
import csv
import glob
import hashlib
import json
import os
import re
import shutil
import sys

(TREE, FROZEN, PLANS_RC8, PLANS_V1, RECORDS, SUITES, AUDITS, TOOLS, RC8_RUNS) = sys.argv[1:10]
sys.path.insert(0, TREE)
from opencif import pipeline as PL                                   # noqa: E402

VAL = os.path.join(TREE, 'validation')
os.makedirs(os.path.join(VAL, 'logs'), exist_ok=True)
frozen = json.load(open(FROZEN))
env = PL.environment_stamp() if hasattr(PL, 'environment_stamp') else PL.ENV

# The rebuild: the last line written for each build is its result.
rebuild = {}
for line in open(os.path.join(TREE, 'runs', 'rebuild_summary.jsonl')):
    d = json.loads(line)
    rebuild[d['name']] = d


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


photomasks, rows, standalone = {}, [], {}
for name in sorted(frozen):
    b = rebuild.get(name)
    if b is None:
        continue
    rec_path = os.path.join(TREE, 'runs', name, name + '.json')
    rec = json.load(open(rec_path)) if os.path.exists(rec_path) else {}
    des, va = rec.get('design') or {}, rec.get('validation') or {}
    pa = va.get('physics_audit') or {}
    sa_path = os.path.join(AUDITS, name + '.json')
    sa = json.load(open(sa_path)) if os.path.exists(sa_path) else None
    photomasks[name] = dict(
        dxf=b.get('dxf'), dxf_sha256=b.get('sha256'), identical_to_rc8=b.get('identical'),
        exit_code=b.get('exit_code'), exit_code_rc8=b.get('exit_code_rc8'), verdict=b.get('verdict'),
        N=des.get('N'), F_design=des.get('F'),
        # a refused build (exit 5) draws no photomask, so neither audit has a file to read
        physics_audit=pa.get('status') or (None if b.get('dxf') else 'not run (no photomask)'),
        physics_audit_standalone=(sa or {}).get('status') or (None if b.get('dxf') else 'not run (no photomask)'),
        code_stamp=b.get('code'),
        record='runs/%s/%s.json' % (name, name))
    for where, a in (('build', pa), ('standalone', sa or {})):
        for c in a.get('checks', []):
            v = c.get('value')
            rows.append(dict(photomask=name, audit=where, check=c['id'], role=c['role'], ok=c['ok'],
                             value=json.dumps(v, sort_keys=True)[:300] if v is not None else ''))

# The arrays, which tools/make_arrays.py composes from built elements, against RC8's arrays.
arrays = {}
for name in sorted(frozen):
    rec_path = os.path.join(TREE, 'runs', name, name + '.json')
    if name in photomasks or not os.path.exists(rec_path):
        continue
    rec = json.load(open(rec_path))
    if 'count' not in rec:
        continue
    rc8 = json.load(open(os.path.join(RC8_RUNS, name, name + '.json')))
    got = sha256(os.path.join(TREE, 'runs', name, rec['dxf'])) if rec.get('dxf') else None
    code = rec.get('code')
    arrays[name] = dict(dxf=rec.get('dxf'), dxf_sha256=got, dxf_sha256_rc8=rc8.get('dxf_sha256'),
                        identical_to_rc8=got is not None and got == rc8.get('dxf_sha256'), count=rec.get('count'),
                        inlet=rec.get('inlet'), element=rec.get('element'),
                        code_stamp=code.get('sha256') if isinstance(code, dict) else code,
                        record='runs/%s/%s.json' % (name, name))

# The plans and the records against RC8's.
p8, p1 = json.load(open(PLANS_RC8)), json.load(open(PLANS_V1))
# Two fields differ by design and are compared apart: the planning time ('seconds'), and the status text of the gate
# P8, which RC8 called a diagnostic and which 1.0.0 audits on every photomask (docs/PHYSICS_AUDIT.md). Every other
# field (each array of the row schedule by sha256, each scalar, the other gates, the options, the notes) must be equal.
PLAN_EXEMPT = ('seconds', 'gates/P8/status')


def plan_core(fp):
    fp = json.loads(json.dumps(fp))
    fp.pop('seconds', None)
    if isinstance(fp.get('gates'), dict) and isinstance(fp['gates'].get('P8'), dict):
        fp['gates']['P8'].pop('status', None)
    return fp


plan_same = sorted(n for n in p8 if n in p1 and plan_core(p8[n]) == plan_core(p1[n]))
plan_diff = sorted(n for n in p8 if n in p1 and plan_core(p8[n]) != plan_core(p1[n]))
plan_missing = sorted(set(p8) ^ set(p1))
plan_p8_status_changed = sum(1 for n in p8 if n in p1 and
                             ((p8[n].get('gates') or {}).get('P8') or {}).get('status') !=
                             ((p1[n].get('gates') or {}).get('P8') or {}).get('status'))
records = json.load(open(RECORDS))

# The suites, from their own logs.
suites = {}
for log in sorted(glob.glob(os.path.join(SUITES, 'test_*.log'))):
    name = os.path.basename(log)[:-4]
    text = open(log, encoding='utf-8', errors='replace').read().splitlines()
    last = text[-1] if text else ''
    m = re.search(r'(\d+) checks', last)
    if m:                                  # a suite that counts its own checks says so on its last line
        n_ok = int(m.group(1))
    else:                                  # otherwise each passed check prints one line marked ok (or a preset OK)
        n_ok = sum(1 for ln in text if re.match(r'^\s*ok\b', ln) or re.match(r'^\S+\s+OK\s', ln)
                   or ' ok ' in ln or ln.rstrip().endswith(' ok'))
    suites[name] = dict(result=last.strip(), checks_ok=n_ok)
    # the tree's absolute folder on the build machine appears in Python's warnings; the copy names it <repo>
    with open(log, encoding='utf-8', errors='replace') as fh:
        body = fh.read().replace(os.path.abspath(TREE) + os.sep, '<repo>/')
    with open(os.path.join(VAL, 'logs', os.path.basename(log)), 'w', encoding='utf-8') as fh:
        fh.write(body)
shutil.copy(os.path.join(TREE, 'runs', 'rebuild_summary.jsonl'), os.path.join(VAL, 'logs', 'rebuild_summary.jsonl'))
# RC8's photomask hashes, both releases' plan fingerprints, and the record comparison, so that every count above can be
# checked from the tree (they hold hashes and numbers, no image)
os.makedirs(os.path.join(VAL, 'rc8_comparison'), exist_ok=True)
for src_, name_ in ((FROZEN, 'frozen_photomasks.json'), (PLANS_RC8, 'plans_rc8.json'), (PLANS_V1, 'plans_1.0.0.json'),
                    (RECORDS, 'compare_records.json')):
    shutil.copy(src_, os.path.join(VAL, 'rc8_comparison', name_))

release = dict(
    version=PL.CODE['version'], code_stamp=PL.CODE['sha256'], code_files=PL.CODE['files'],
    tests_stamp=env['tests_sha256'], python=env['python'], packages=env['installed'],
    rc8_code_stamp='593d977fc83ec217',
    counts=dict(photomasks=len(photomasks),
                builds=len(photomasks), drawn=sum(1 for v in photomasks.values() if v['dxf']),
                refused=sum(1 for v in photomasks.values() if v['exit_code'] == 5),
                arrays=len(arrays), arrays_identical_to_rc8=sum(1 for v in arrays.values() if v['identical_to_rc8']),
                identical_to_rc8=sum(1 for v in photomasks.values() if v['identical_to_rc8']),
                exit_codes_as_rc8=sum(1 for v in photomasks.values() if v['exit_code'] == v['exit_code_rc8']),
                audit_in_build=dict((s, sum(1 for v in photomasks.values() if v['physics_audit'] == s))
                                    for s in sorted({v['physics_audit'] for v in photomasks.values()}, key=str)),
                audit_standalone=dict((s, sum(1 for v in photomasks.values() if v['physics_audit_standalone'] == s))
                                      for s in sorted({v['physics_audit_standalone'] for v in photomasks.values()},
                                                      key=str)),
                plans_compared=len(plan_same) + len(plan_diff), plans_identical=len(plan_same),
                plans_different=plan_diff, plans_missing=plan_missing, plans_exempt_fields=list(PLAN_EXEMPT),
                plans_p8_status_changed=plan_p8_status_changed, records=records.get('summary', records) if isinstance(records, dict)
                else records),
    suites=suites, photomasks=photomasks, arrays=arrays)
json.dump(release, open(os.path.join(VAL, 'RELEASE.json'), 'w'), indent=1, sort_keys=True, default=str)
with open(os.path.join(VAL, 'audit_table.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['photomask', 'audit', 'check', 'role', 'ok', 'value'])
    w.writeheader()
    w.writerows(rows)
with open(os.path.join(VAL, 'photomasks_sha256.txt'), 'w') as fh:
    for name, v in sorted(list(photomasks.items()) + list(arrays.items())):
        if v['dxf']:
            fh.write('%s  %s/%s\n' % (v['dxf_sha256'], name, v['dxf']))
# The restructure: every script that turns RC8's frozen tree into this one, with the voiced modules and documents.
rs = os.path.join(VAL, 'restructure')
if os.path.exists(rs):
    shutil.rmtree(rs)
for root, dirs, files in os.walk(TOOLS):           # every file that make_tree.sh reads (no caches; prevoice/ holds the
    dirs[:] = [d for d in dirs if d not in ('__pycache__', 'prevoice')]     # copies before the voice pass, not read)
    for fn in files:
        if not fn.endswith('.pyc'):
            src = os.path.join(root, fn)
            dst = os.path.join(rs, os.path.relpath(src, TOOLS))
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy(src, dst)
print('validation/: %d builds (%d identical to RC8), %d arrays (%d identical), %d audit rows, %d suites, %d plans compared'
      % (len(photomasks), release['counts']['identical_to_rc8'], len(arrays),
         release['counts']['arrays_identical_to_rc8'], len(rows), len(suites), release['counts']['plans_compared']))
