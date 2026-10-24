"""rc7 (code review R12): a fingerprint of the planner and the head geometry of a spec -- the schedule (S.plan: kinds,
x, ell, w_c, f, w_s) and every head unit as drawn (geometry.plan_head, solve_transition, exact_unit, unit_path's
points) -- so that a test can show that rc7 leaves rc6's arc-turn heads unchanged without importing rc6.  The same
script runs on either tree (its ROOT argument).
usage: python3 -I fingerprint_head.py ROOT SPEC.yaml [SPEC.yaml ...]   -> one JSON line {spec: sha256[:16]}"""
import hashlib
import json
import os
import sys

ROOT = os.path.abspath(sys.argv[1])
sys.path[:0] = [ROOT]
from opencif import spec as SP                              # noqa: E402
from opencif.physics import schedule as S  # noqa: E402
from opencif.geometry import shapes as Gm  # noqa: E402

out = {}
for f in sys.argv[2:]:
    sp, _ = SP.to_engine(SP.load(f if os.path.isabs(f) else os.path.join(ROOT, f)))
    S.clear_caches() if hasattr(S, 'clear_caches') else (S._BEND_CACHE.clear(), S._LAM_CACHE.clear())
    sc = S.plan(sp)
    rec = dict(kind=[str(k) for k in sc.kind], x=[float(v) for v in sc.x], ell=[float(v) for v in sc.ell],
               wc=[float(v) for v in sc.wc], f=[float(v) for v in sc.f], ws=[float(v) for v in sc.ws])
    units = Gm.plan_head(sc)
    if units:
        Gm.solve_transition(sc, units[-1])
    for u in units:
        Gm.exact_unit(sp, u)
    rec['units'] = [[[float(a), float(b)] for a, b in Gm.unit_path(sp, u).p] for u in units]
    h = hashlib.sha256(json.dumps(rec, sort_keys=True).encode()).hexdigest()[:16]
    out[os.path.basename(f)] = h
print(json.dumps(out, sort_keys=True))
