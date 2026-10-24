"""The Dinh V1 presets' OWN schedule (the printed 75 -> 52 um centre channel, the assumed constant-shear taper, the
reconstructed head, the printed f*'s lamina held at every gap) planned to other filtrate fractions F, with the length
it needs: the head and the linear rows (and, given the presets' run records, the drawn parts outside them: inlet,
filter, outlet).  The schedule-specific companion of dinh_length_bound.py (which bounds every taper).

rc6 port of provenance/rc5_scripts/dinh_schedule_at_F.py: the head's upstream end x_top comes from the plan itself
(layout.head_top_x, as the build computes it), so the head and rows need no build; with --runs the record's drawn
footprint adds the other parts, and its x_top is checked against the plan's.  The row pitch is the preset's (rc6
45.1 um).  On the rc5 tree with --runs it reproduces rc5's dinh_schedule_at_F.json
(provenance/rc6_scripts/out/length_regression_rc5.out).
Usage: python3 dinh_schedule_at_F.py OPENCIF_ROOT OUT.json [--runs RUNS_DIR] [F ...]
       (default F 0.5592 0.5105 0.505 0.472)"""
import copy
import json
import os
import sys

argv, RUNS = sys.argv[1:], None
if '--runs' in argv:
    i_ = argv.index('--runs')
    RUNS = argv[i_ + 1]
    del argv[i_:i_ + 2]
R, OUT = argv[:2]
Fs = [float(v) for v in argv[2:]] or [0.5592, 0.5105, 0.505, 0.472]
sys.path[:0] = [R]
from opencif import spec as SP
from opencif.physics import schedule as S
from opencif.geometry import layout as L

out = {}
for name in ('dinh2024_v1_depicted', 'dinh2024_v1_specified'):
    rec = None
    if RUNS is not None:
        rec = json.load(open(os.path.join(RUNS, name, name + '.json')))
    entry = dict(plans=[])
    if rec is not None:
        dp = rec['validation']['drawn_parameters']
        x_top_rec = rec['design']['x_top']
        other_mm = dp['footprint_mm'][0] - (dp['array_length_mm'] - x_top_rec / 1000.0)   # inlet, filter, outlet parts
        entry.update(drawn_footprint_mm=dp['footprint_mm'][0], drawn_other_parts_mm=round(other_mm, 3))
    out[name] = entry
    for F in Fs:
        d = copy.deepcopy(SP.load(os.path.join(R, 'presets', name + '.yaml')))
        d['targets']['filtrate_fraction'] = F
        sp, op = SP.to_engine(SP.load(d))
        sc = S.plan(sp)
        N, NM = sc.N, sc.NM
        x_top = float(L.head_top_x(sc))
        if rec is not None:
            assert abs(x_top - x_top_rec) < 1e-6, (name, F, x_top, x_top_rec)
        head_mm = (sc.x[NM] - x_top) / 1000.0 if NM else 0.0
        lin_mm = (sc.x[N - 1] - sc.x[NM]) / 1000.0
        row = dict(F_target=F, F=round(float(sc.F), 5), N=int(N), NM=int(NM), linear_rows=int(N - NM),
                   wc_last_um=round(float(sc.wc[-1]), 3), head_mm=round(head_mm, 3), linear_rows_mm=round(lin_mm, 3),
                   head_plus_rows_mm=round(head_mm + lin_mm, 3))
        if rec is not None:
            row['with_other_parts_mm'] = round(head_mm + lin_mm + other_mm, 3)
        out[name]['plans'].append(row)
        print(name, json.dumps(row), flush=True)
json.dump(out, open(OUT, 'w'), indent=1)
