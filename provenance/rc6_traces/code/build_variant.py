"""Build a variant spec (not a preset) into its own directory: python3 build_variant.py ROOT SPEC.yaml OUT_DIR
(pipeline.run with that output directory: the DXF, record, pickle and figures there, nothing in runs/)."""
import sys

R, SPEC, OUT = sys.argv[1:4]
sys.path[:0] = [R]
from opencif import pipeline as PL
rec = PL.run(SPEC, outdir=OUT, verbose=True, render=True)
a = rec.get('assessment') or {}
print('exit', a.get('exit_code'), a.get('verdict', '')[:300])
sys.exit(0 if a.get('exit_code') in (0, 3) else 1)
