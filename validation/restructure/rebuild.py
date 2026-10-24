"""Rebuild every frozen photomask with a given tree, in lanes, resumably, within a deadline.

usage: python -I rebuild.py TREE FROZEN_JSON RC8_RUNS LANES DEADLINE_EPOCH [NAME ...]
  TREE          the opencif tree to build with (builds go to TREE/runs/<name>/, as tools/ expects)
  FROZEN_JSON   rc8ref/frozen_photomasks.json (name -> the RC8 photomask's sha256, exit code)
  RC8_RUNS      RC8's runs/ (for each build's duration, to schedule longest first and within the deadline)
  LANES         parallel builds
  DEADLINE      wall-clock epoch seconds: a build is started only if it is expected to finish before it
A build whose folder already holds a complete record of this tree's code stamp is skipped (resumable).  One JSON line
per finished build goes to TREE/runs/rebuild_summary.jsonl."""
import glob
import hashlib
import json
import os
import subprocess
import sys
import time

TREE, FROZEN, RC8, LANES, DEADLINE = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5])
ONLY = set(sys.argv[6:])
PY = sys.executable
frozen = json.load(open(FROZEN))
specs = {os.path.splitext(os.path.basename(p))[0]: p for p in
         glob.glob(os.path.join(TREE, 'presets', '*.yaml')) + glob.glob(os.path.join(TREE, 'examples', '*.yaml')) +
         glob.glob(os.path.join(TREE, 'examples', '*', '*.yaml')) +
         glob.glob(os.path.join(TREE, 'provenance', 'rc7_traces', 'variants', '*.yaml'))}
stamp = subprocess.run([PY, '-I', '-B', '-c', 'import sys; sys.path.insert(0, %r); from opencif import pipeline as P; '
                        'print(P.CODE["sha256"])' % TREE], capture_output=True, text=True).stdout.strip()
RUNS = os.path.join(TREE, 'runs')
os.makedirs(RUNS, exist_ok=True)
SUMMARY = os.path.join(RUNS, 'rebuild_summary.jsonl')


def expected(name):
    try:
        r = json.load(open(os.path.join(RC8, name, name + '.json')))
        return float(r.get('seconds_total') or r.get('seconds') or 600.0)
    except Exception:
        return 600.0


def complete(name):
    p = os.path.join(RUNS, name, name + '.json')
    try:
        r = json.load(open(p))
    except Exception:
        return False
    return (r.get('code') or {}).get('sha256') == stamp and ('files' in r or (r.get('assessment') or {}).get('exit_code') == 5)


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest() if os.path.exists(p) else None


todo = [n for n in frozen if n in specs and (not ONLY or n in ONLY) and not complete(n)]
todo.sort(key=expected, reverse=True)                                   # longest first
print('code %s: %d to build (%d already complete); deadline in %.0f min' % (
    stamp, len(todo), sum(1 for n in frozen if n in specs and complete(n)), (DEADLINE - time.time()) / 60), flush=True)
running = {}
while todo or running:
    now = time.time()
    while len(running) < LANES:
        nxt = next((n for n in todo if now + 1.25 * expected(n) + 60 < DEADLINE), None)
        if nxt is None:
            break
        todo.remove(nxt)
        od = os.path.join(RUNS, nxt)
        os.makedirs(od, exist_ok=True)
        lg = open(os.path.join(od, 'build.log'), 'w')
        cmd = [PY, '-I', '-B', '-c', 'import sys; sys.path.insert(0, %r); from opencif.__main__ import main; '
               'sys.exit(main(["build", %r, "--out", %r]))' % (TREE, specs[nxt], od)]
        running[nxt] = (subprocess.Popen(cmd, cwd=TREE, stdout=lg, stderr=subprocess.STDOUT), time.time(), lg)
        print('%s start %s (expected %.0f s)' % (time.strftime('%H:%M:%S'), nxt, expected(nxt)), flush=True)
    if not running:
        break                                                            # nothing fits before the deadline
    time.sleep(5)
    for name, (p, ts, lg) in list(running.items()):
        if p.poll() is None:
            continue
        lg.close()
        del running[name]
        od = os.path.join(RUNS, name)
        try:
            r = json.load(open(os.path.join(od, name + '.json')))
        except Exception:
            r = {}
        a = r.get('assessment') or {}
        files = r.get('files') or {}
        dxf = os.path.join(od, files['dxf']) if files.get('dxf') else None
        got = sha(dxf) if dxf else None
        fz = frozen[name]
        pa = (r.get('validation') or {}).get('physics_audit') or {}
        line = dict(name=name, rc=p.returncode, seconds=round(time.time() - ts), exit_code=a.get('exit_code'),
                    exit_code_rc8=fz.get('exit_code'), dxf=files.get('dxf'), dxf_rc8=fz.get('dxf'), sha256=got,
                    identical=(got == fz.get('dxf_sha256')) if fz.get('dxf_sha256') else (got is None),
                    audit=pa.get('status'), audit_failed=pa.get('hard_failures'), audit_envelope=pa.get('envelope_failures'),
                    code=(r.get('code') or {}).get('sha256'), verdict=(a.get('verdict') or '')[:200])
        with open(SUMMARY, 'a') as fh:
            fh.write(json.dumps(line) + '\n')
        print('%s done  %s rc %s exit %s (rc8 %s) identical %s audit %s %s (%d s)' % (
            time.strftime('%H:%M:%S'), name, p.returncode, line['exit_code'], line['exit_code_rc8'], line['identical'],
            line['audit'], line['audit_failed'] or '', line['seconds']), flush=True)
print('left for the next window: %s' % todo, flush=True)
