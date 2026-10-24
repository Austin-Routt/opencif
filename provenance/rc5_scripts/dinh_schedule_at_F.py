"""rc5 (after the rc5 reviews): the Dinh V1 presets' OWN schedule (the printed 75 -> 52 um centre channel, the assumed
constant-shear taper, the reconstructed head, the printed f*'s lamina held at every gap) planned to other filtrate fractions F,
with the length it needs: the head, the linear rows, and the drawn preset's other parts (inlet, filter, outlet) carried
over from its record.  The schedule-specific companion of dinh_length_bound.py (which bounds every taper).
Usage: python3 dinh_schedule_at_F.py OPENCIF_ROOT [F ...]   (default F 0.5592 0.5105 0.505 0.472)"""
import copy, json, os, sys
R = sys.argv[1]
Fs = [float(v) for v in sys.argv[2:]] or [0.5592, 0.5105, 0.505, 0.472]
sys.path[:0] = [R, os.path.join(R, 'opencif'), os.path.join(R, 'opencif', 'engine')]
from opencif import spec as SP
import cifsched as S
out = {}
for name in ('dinh2024_v1_depicted', 'dinh2024_v1_specified'):
    rec = json.load(open(os.path.join(R, 'runs', name, name + '.json')))
    dp = rec['validation']['drawn_parameters']
    x_top = rec['design']['x_top']
    other_mm = dp['footprint_mm'][0] - (dp['array_length_mm'] - x_top / 1000.0)   # inlet, filter, outlet parts
    out[name] = dict(drawn_footprint_mm=dp['footprint_mm'][0], drawn_other_parts_mm=round(other_mm, 3), plans=[])
    for F in Fs:
        d = copy.deepcopy(SP.load(os.path.join(R, 'presets', name + '.yaml')))
        d['targets']['filtrate_fraction'] = F
        sp, op = SP.to_engine(SP.load(d))
        sc = S.plan(sp)
        N, NM = sc.N, sc.NM
        head_mm = (sc.x[NM] - x_top) / 1000.0 if NM else 0.0
        lin_mm = (sc.x[N - 1] - sc.x[NM]) / 1000.0
        row = dict(F_target=F, F=round(float(sc.F), 5), N=int(N), NM=int(NM), linear_rows=int(N - NM),
                   wc_last_um=round(float(sc.wc[-1]), 3), head_mm=round(head_mm, 3), linear_rows_mm=round(lin_mm, 3),
                   head_plus_rows_mm=round(head_mm + lin_mm, 3), with_other_parts_mm=round(head_mm + lin_mm + other_mm, 3))
        out[name]['plans'].append(row)
        print(name, json.dumps(row), flush=True)
json.dump(out, open('dinh_schedule_at_F.json', 'w'), indent=1)
