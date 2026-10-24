"""Plan every shipped specification (presets, examples, the Fig. 4 variants) in one tree and write a fingerprint of
everything the plan computes: every array of the row schedule (sha256 of its bytes) and every scalar (repr), and the
gates.  Run it in the RC8 tree and in the 1.0.0 tree and compare: the restructure must change nothing.
usage: python -I regress_plans.py TREE OUT.json [GLOB_FILTER]"""
import glob
import hashlib
import json
import os
import sys
import time

TREE, OUT = sys.argv[1:3]
FILT = sys.argv[3] if len(sys.argv) > 3 else ''
sys.path.insert(0, TREE)
from opencif import pipeline as PL, spec as SP   # noqa: E402
try:
    from opencif.physics import schedule as S    # 1.0.0
except ImportError:
    import cifsched as S                          # RC8 (the engine on sys.path)
import numpy as np                                # noqa: E402


def fp(v, depth=0):
    if isinstance(v, np.ndarray):
        return 'nd:%s:%s:%s' % (v.dtype, v.shape, hashlib.sha256(np.ascontiguousarray(v).tobytes()).hexdigest()[:24])
    if isinstance(v, (list, tuple)):
        if v and all(isinstance(x, (int, float, np.floating, np.integer)) for x in v):
            return 'seq:' + hashlib.sha256(repr([float(x) for x in v]).encode()).hexdigest()[:24]
        return [fp(x, depth + 1) for x in v] if depth < 3 else repr(v)[:200]
    if isinstance(v, dict):
        return {str(k): fp(x, depth + 1) for k, x in v.items()} if depth < 4 else repr(v)[:200]
    if isinstance(v, (float, int, str, bool, type(None), np.floating, np.integer)):
        return repr(v)
    if hasattr(v, '__dict__') and depth < 3:
        return {k: fp(x, depth + 1) for k, x in vars(v).items() if not k.startswith('_')}
    return type(v).__name__


specs = sorted(glob.glob(os.path.join(TREE, 'presets', '*.yaml')) + glob.glob(os.path.join(TREE, 'examples', '*.yaml')) +
               glob.glob(os.path.join(TREE, 'examples', '*', '*.yaml')) +
               glob.glob(os.path.join(TREE, 'provenance', 'rc7_traces', 'variants', '*.yaml')))
res = {}
for p in specs:
    rel = os.path.relpath(p, TREE)
    if FILT and FILT not in rel:
        continue
    t0 = time.time()
    try:
        S.clear_caches()
        spec, opts, sc, g, env_bad, notes = PL.prepare(SP.load(p), lambda m: None)
        res[rel] = dict(spec=fp(spec), sc=fp(sc), gates=fp(g), env_bad=list(env_bad), notes=[str(n) for n in notes],
                        opts=fp(opts), seconds=round(time.time() - t0, 1))
    except Exception as e:          # a refused or invalid spec is a result too
        res[rel] = dict(error='%s: %s' % (type(e).__name__, str(e)[:300]), seconds=round(time.time() - t0, 1))
    print('%-75s %6.1f s %s' % (rel, time.time() - t0, 'ERROR' if 'error' in res[rel] else ''), flush=True)
json.dump(res, open(OUT, 'w'), indent=1, sort_keys=True)
print('planned', len(res))
