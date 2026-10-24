"""rc7 (code review R2): regenerate data/depicted.json['dinh'] -- the gap fractions that make the closure draw Dinh
Fig. 2A inset's units -- for the head as it is now drawn: square quarter-turns with the measured roundings (outer 7 um,
inner 10 um; DINH-1), instead of rc6's arc turns.  Every other field of the spec is devspecs.dinh_fig2('V1'), which equals
the engine spec of presets/dinh2024_v1_depicted.yaml in every field the head uses (checked below).  data/depicted.json
['xia'] is copied unchanged (this script does not recompute it).
usage: python3 -I provenance/rc7_scripts/regen_depicted_dinh_square.py [--check]   (--check: compare, write nothing)"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from opencif import spec as SP                              # noqa: E402  (sets the engine path)
from opencif import published as D  # noqa: E402
from opencif.geometry import depicted as DP  # noqa: E402

OUT = os.path.join(ROOT, 'opencif', 'data', 'depicted.json')
sp = D.dinh_fig2('V1')
sp = sp.__class__(**{**sp.__dict__, 'turns': 'square', 'r_corner_out': 7.0, 'r_corner_in': 10.0})
# the preset's engine spec, for the fields the head uses
pre, _ = SP.to_engine(SP.load(os.path.join(ROOT, 'presets', 'dinh2024_v1_depicted.yaml')))
skip = {'f_head', 'name', 'note', 'flare_deg', 'r_junction', 'exp_domain_refuse', 'tubing_centre_id_mm'}
diff = sorted(k for k in set(sp.__dict__) | set(pre.__dict__) if k not in skip and sp.__dict__.get(k) != pre.__dict__.get(k))
assert not diff, ('the preset differs from the regeneration spec in', diff)
res = DP.report('Dinh 2024 Fig. 2A inset (V1), square quarter-turns (outer 7 um, inner 10 um; rc7 DINH-1, code review R2)',
                sp, DP.dinh_depicted_geometry(sp))
res['turns'] = dict(turns='square', r_corner_out_um=7.0, r_corner_in_um=10.0,
                    note='rc7 (code review R2): regenerated for the square-turn head (rc6 derived it for arc turns; '
                         'see provenance/rc7_scripts/regen_depicted_dinh_square.py)')
with open(OUT) as fh:
    old = json.load(fh)
for a, b in zip(old['dinh']['f_head'], res['f_head']):
    print('f %.6e -> %.6e  (%+.3f %%)' % (a, b, 100.0 * (b / a - 1.0)))
if '--check' not in sys.argv:
    new = dict(old)
    new['dinh'] = res
    with open(OUT, 'w') as fh:
        json.dump(new, fh, indent=1)
    with open(OUT) as fh:
        chk = json.load(fh)
    assert json.dumps(chk['xia'], sort_keys=True) == json.dumps(old['xia'], sort_keys=True)
    print('written', OUT, '(xia unchanged)')
