"""Tally every changed leaf path between the rc3 and rc4 records (lists of dicts keyed by 'check' where present)."""
import glob, json, os, collections, sys
A = sys.argv[1] if len(sys.argv) > 2 else 'opencif_next/runs'   # rc3 runs/
B = sys.argv[2] if len(sys.argv) > 2 else 'opencif_rc4/runs'    # rc4 runs/
def load(root, n):
    p = os.path.join(root, n, n + '.json'); return json.load(open(p)) if os.path.exists(p) else None
def norm(o):
    if isinstance(o, list) and o and all(isinstance(x, dict) and 'check' in x for x in o):
        return {('CHK:' + x['check'][:38]): x for x in o}
    return o
def walk(a, b, path, out):
    a, b = norm(a), norm(b)
    if isinstance(a, dict) and isinstance(b, dict):
        for k in set(a) | set(b):
            if k not in a: out.append((path + '/' + str(k), 'ADDED'))
            elif k not in b: out.append((path + '/' + str(k), 'REMOVED'))
            else: walk(a[k], b[k], path + '/' + str(k), out)
    elif a != b:
        if isinstance(a, float) and isinstance(b, float) and abs(a - b) <= 1e-12 * max(1, abs(a)): return
        out.append((path, 'CHANGED'))
names = sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(A, '*', '*.json'))
               if os.path.basename(os.path.dirname(p)) + '.json' == os.path.basename(p))
tally = collections.defaultdict(list)
for n in names:
    a, b = load(A, n), load(B, n); out = []
    walk(a, b, '', out)
    for p, kind in out:
        if '/code/' in p + '/' or p.startswith('/files/') or p.startswith('/created') or 'element_code' in p: continue
        tally[(p, kind)].append(n)
for (p, kind), ns in sorted(tally.items(), key=lambda t: (-len(t[1]), t[0])):
    print('%3d %-8s %s   %s' % (len(ns), kind, p[:150], (', '.join(ns) if len(ns) <= 8 else '')))
