"""run_batch.py: build a list of specs in LANES parallel lanes.

Usage: python3 tools/run_batch.py LANES spec1.yaml spec2.yaml ...

The specs start in the order given, because the script does not sort them, so the caller lists the longest first. Each
build writes its log to logs/<name>.log, and logs/batch_summary.txt is rewritten as each build ends."""
import os, subprocess, sys, time
lanes = int(sys.argv[1]); specs = sys.argv[2:]
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
queue = list(specs); running = {}; done = []
t0 = time.time()
while queue or running:
    while queue and len(running) < lanes:
        sp = queue.pop(0); name = os.path.splitext(os.path.basename(sp))[0]
        lg = open(os.path.join(ROOT, 'logs', name + '.log'), 'w')
        running[name] = (subprocess.Popen([sys.executable, '-m', 'opencif', 'build', sp], cwd=ROOT, stdout=lg, stderr=subprocess.STDOUT), time.time(), lg)
    time.sleep(5)
    for name, (p, ts, lg) in list(running.items()):
        if p.poll() is not None:
            lg.close()
            tail = open(os.path.join(ROOT, 'logs', name + '.log')).read().strip().splitlines()[-1:]
            done.append((name, p.returncode, round(time.time() - ts), tail[0][:300] if tail else ''))
            del running[name]
            with open(os.path.join(ROOT, 'logs', 'batch_summary.txt'), 'w') as fh:
                for d in done:
                    fh.write('%-44s rc %s %5ds  %s\n' % d)
print('batch done in %.0f s' % (time.time() - t0))
