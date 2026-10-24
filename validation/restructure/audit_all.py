"""Audit every frozen photomask from scratch with a given tree (opencif.checks.auditor.audit_run: plan and draw again,
read the DXF back, measure the file's polygon), in lanes, resumably, within a deadline.
usage: python -I audit_all.py TREE RUNS_DIR FROZEN_JSON OUT_DIR LANES DEADLINE_EPOCH [--drc] [NAME ...]"""
import json
import os
import subprocess
import sys
import time

args = [a for a in sys.argv[1:] if a != '--drc']
DRC = '--drc' in sys.argv
TREE, RUNS, FROZEN, OUT, LANES, DEADLINE = args[0], args[1], args[2], args[3], int(args[4]), float(args[5])
ONLY = set(args[6:])
PY = sys.executable
os.makedirs(OUT, exist_ok=True)
frozen = json.load(open(FROZEN))


def expected(n):
    try:
        r = json.load(open(os.path.join(RUNS, n, n + '.json')))
        return float(r.get('seconds_total') or r.get('seconds') or 600.0) * 0.6
    except Exception:
        return 400.0


names = [n for n in frozen if os.path.exists(os.path.join(RUNS, n, n + '.json')) and (not ONLY or n in ONLY)
         and not os.path.exists(os.path.join(OUT, n + '.json'))]
names = [n for n in names if (json.load(open(os.path.join(RUNS, n, n + '.json'))).get('files') or {}).get('dxf')]
names.sort(key=expected, reverse=True)
print('%d to audit; deadline in %.0f min' % (len(names), (DEADLINE - time.time()) / 60), flush=True)
prog = ('import sys, json; sys.path.insert(0, %r); from opencif.checks import auditor as AU\n'
        'a = AU.audit_run(sys.argv[1], sys.argv[2], say=lambda m: print(m, flush=True), drc=%r)\n'
        'json.dump(a, open(sys.argv[3], "w"), indent=1, default=float)\n') % (TREE, DRC)
running = {}
while names or running:
    now = time.time()
    while len(running) < LANES:
        nxt = next((n for n in names if now + 1.5 * expected(n) + 60 < DEADLINE), None)
        if nxt is None:
            break
        names.remove(nxt)
        lg = open(os.path.join(OUT, nxt + '.log'), 'w')
        running[nxt] = (subprocess.Popen([PY, '-I', '-B', '-c', prog, os.path.join(RUNS, nxt), nxt,
                                          os.path.join(OUT, nxt + '.json.part')], stdout=lg, stderr=subprocess.STDOUT),
                        time.time(), lg)
        print('%s start %s' % (time.strftime('%H:%M:%S'), nxt), flush=True)
    if not running:
        break
    time.sleep(5)
    for n, (p, ts, lg) in list(running.items()):
        if p.poll() is None:
            continue
        lg.close()
        del running[n]
        part = os.path.join(OUT, n + '.json.part')
        line = dict(name=n, rc=p.returncode, seconds=round(time.time() - ts))
        if p.returncode == 0 and os.path.exists(part):
            os.replace(part, os.path.join(OUT, n + '.json'))
            a = json.load(open(os.path.join(OUT, n + '.json')))
            ck = {c['id']: c for c in a['checks']}
            line.update(status=a['status'], hard=a['hard_failures'], envelope=a['envelope_failures'],
                        R0=ck['R0']['value'].get('max_rel_diff'), R0_segments=ck['R0']['value'].get('segments'),
                        P4=ck['P4']['value'].get('bulk_max'), P5=ck['P5']['value'].get('max_delta_over_delta_star'),
                        R6b=ck['R6b']['value'].get('max_dev_across_um'))
        with open(os.path.join(OUT, 'summary.jsonl'), 'a') as fh:
            fh.write(json.dumps(line) + '\n')
        print('%s done  %s' % (time.strftime('%H:%M:%S'), json.dumps(line)), flush=True)
print('left: %s' % names, flush=True)
