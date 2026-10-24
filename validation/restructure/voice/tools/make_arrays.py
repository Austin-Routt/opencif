"""make_arrays.py: compose a parallel array from a finished, validated element run, without rebuilding the element.

Usage: python3 tools/make_arrays.py RUN_NAME COUNT [shared|separate] [chip_width_mm chip_length_mm]

The array is written to runs/<RUN_NAME>_x<COUNT>/."""
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from opencif import compose as CP
name, count = sys.argv[1], int(sys.argv[2])
inlet = sys.argv[3] if len(sys.argv) > 3 else 'shared'
chip = (float(sys.argv[5]) * 1000.0, float(sys.argv[4]) * 1000.0) if len(sys.argv) > 5 else None
d = os.path.join(ROOT, 'runs', name)
od = os.path.join(ROOT, 'runs', '%s_x%d' % (name, count))
# This array's earlier outputs are removed first, as pipeline.run does for a build, because pipeline.name_mask does not
# replace an existing <name>.dxf, so an array composed into a used folder would keep the earlier photomask under the
# deliverable name (LIMITATIONS section E; final document check, M1).
if os.path.isdir(od):
    import re
    own = re.compile(re.escape('%s_x%d' % (name, count)) + r'(\.(json|dxf|png)|_[a-z]+\.(json|png)|_FAILED\.dxf|_UNASSESSED\.dxf)$')
    for f in os.listdir(od):
        if own.match(f):
            os.remove(os.path.join(od, f))
out = CP.array(d, name, count, od, inlet=inlet, chip=chip)
import json
print(json.dumps({k: v for k, v in out.items() if k != 'checks'}, indent=1, default=float))
print(json.dumps(out['checks'], indent=1, default=float))
