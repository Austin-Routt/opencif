"""rc5: the final rebuild against the first full rc5 rebuild (code f6c8eba91ec1f316), record by record and mask by mask.
Tallies every changed JSON field (timings, logs and the code stamp aside) and compares each DXF's sha256.
Usage: python3 rebuild_diff.py OLD_RUNS NEW_RUNS   (OLD_RUNS may hold <name>/dxf.sha256 instead of the DXFs)"""
import collections, glob, hashlib, json, os, re, sys
A, B = sys.argv[1], sys.argv[2]
SKIP = re.compile(r'(seconds|/log$|^/code|/element_code|/files/|/timing)')


def walk(a, b, path=''):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            out += walk(a.get(k), b.get(k), path + '/' + str(k))
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b) and all(not isinstance(x, (dict, list)) for x in a + b) is False:
        for i, (x, y) in enumerate(zip(a, b)):
            out += walk(x, y, path + '[]')
    elif a != b:
        out.append(path)
    return out


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


tally, per, masks = collections.Counter(), collections.defaultdict(list), {'same': [], 'changed': [], 'missing': []}
for jp in sorted(glob.glob(os.path.join(B, '*', '*.json'))):
    n = os.path.basename(os.path.dirname(jp))
    if os.path.basename(jp) != n + '.json':
        continue
    op = os.path.join(A, n, n + '.json')
    if not os.path.exists(op):
        tally['NO OLD RECORD'] += 1
        continue
    for p in set(walk(json.load(open(op)), json.load(open(jp)))):
        if SKIP.search(p):
            continue
        tally[p] += 1
        per[p].append(n)
    for dx in glob.glob(os.path.join(B, n, '*.dxf')):
        old = None
        if os.path.exists(os.path.join(A, n, os.path.basename(dx))):
            old = sha(os.path.join(A, n, os.path.basename(dx)))
        elif os.path.exists(os.path.join(A, n, 'dxf.sha256')):
            for line in open(os.path.join(A, n, 'dxf.sha256')):
                h, f = line.split()
                if f == os.path.basename(dx):
                    old = h
        (masks['missing'] if old is None else masks['same'] if old == sha(dx) else masks['changed']).append(n)
for p, c in sorted(tally.items(), key=lambda x: (-x[1], x[0])):
    print('%3d  %-90s %s' % (c, p[:90], ', '.join(per[p][:4]) + (' ...' if len(per[p]) > 4 else '')))
win = os.path.join(A, 'cfd_windows.sha256')                  # the CFD windows: <hash>  runs/<name>/cfd_kit_*/window.dxf
if os.path.exists(win):
    for line in open(win):
        h, f = line.split()
        new = os.path.join(B, os.path.relpath(f, 'runs'))
        (masks['missing'] if not os.path.exists(new) else masks['same'] if sha(new) == h else masks['changed']).append(f)
print('MASKS (designs, arrays, CFD windows) same %d, changed %d %s, missing %d %s' % (
    len(masks['same']), len(masks['changed']), masks['changed'], len(masks['missing']), masks['missing']))
