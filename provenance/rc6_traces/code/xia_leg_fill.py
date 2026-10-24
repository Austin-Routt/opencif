"""rc6 (the rc6 traces report, sec. 12.4): the Xia RS-7 length statement.  Fig. 3a's legs, read at Fig. 3a's own 5 mm
bar (178.32 px, its anti-aliased ends), against the rows the presets' schedules place in them at the row pitch measured
on Fig. 3's panels c and g (51.78 um at 1.798 um/px), and what the model would need to fill them.

Model-free part: the straight leg between loop centres on Fig. 3a (opencif/data/r4/meas_xia_fig3a.json, read at
179 px, x the rc6 factor 178.32 -> 179) and on each rc6 record (validation.length_budget.leg_path_lengths_mm, legs 2-6,
whose rows run from loop centre to loop centre: rows x pitch); Fig. 3a's extent against the drawn footprint.

Conditional part (the model's premises: the frozen cutoff law's cubic kernel with the lamina the same at every gap,
the plateau prefactor, D fixed by the printed f_ref through corrected eqn 1, F by the printed ~15:1, the 300 -> 150 um
constant-shear taper): the schedule re-planned (cifsched.plan, the presets' own spec) at other F with D fixed, and at
other D with F fixed, until its rows grow by Fig. 3a's ratio of full-leg rows to the record's.  That ratio is applied
to every leg alike (an approximation: leg 1 holds the head and leg 7 reaches the ring on both).

Usage: python3 xia_leg_fill.py OPENCIF_ROOT RUNS_DIR OUT.json [preset ...]
       (default both: xia2016_rs7_specified xia2016_rs7_depicted; each plan takes 1-2 min)"""
import copy
import json
import math
import os
import sys

R, RUNS, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
NAMES = sys.argv[4:] or ['xia2016_rs7_specified', 'xia2016_rs7_depicted']
sys.path[:0] = [R]
from opencif import spec as SP
from opencif.physics import schedule as S

m6 = json.load(open(os.path.join(R, 'provenance', 'rc6_scripts', 'rc6_measurements.json')))['xia']
m3a = json.load(open(os.path.join(R, 'opencif', 'data', 'r4', 'meas_xia_fig3a.json')))
f3a = m6['fig3a_length_factor_vs_rc5']
Lam = m6['row_pitch_um']
fig = dict(bar_note='Fig. 3a at its 5 mm bar read with its anti-aliased ends (rc4 measurement x %.6f)' % f3a,
           row_pitch_um=Lam, row_pitch_basis='Fig. 3 panels c and g at %.3f um/px' % m6['scale_used_um_per_px'],
           leg_straight_um=m6['fig3a_leg_straight_length_um'],
           rows_per_full_leg=m6['fig3a_leg_straight_length_um'] / Lam,
           extent_mm=[m3a['extent_um']['length'] * f3a / 1000.0, m3a['extent_um']['width'] * f3a / 1000.0])
# the same legs at the panel-i bar's scale (the scale rc5 used for panels b-g): for comparison only
Lam_bar = Lam * m6['scale_um_per_px']['bar_panel_i'] / m6['scale_used_um_per_px']
fig['at_panel_i_bar_scale'] = dict(row_pitch_um=Lam_bar, rows_per_full_leg=m6['fig3a_leg_straight_length_um'] / Lam_bar)
out = dict(fig3a=fig, presets={})
print(json.dumps(fig), flush=True)


def plan(d, F=None, D=None):
    d = copy.deepcopy(d)
    if F is not None:
        d['targets']['filtrate_fraction'] = F
    if D is not None:
        d['targets']['cutoff_um'] = D
    sp, _ = SP.to_engine(SP.load(d))
    sc = S.plan(sp)
    return int(sc.N), float(sc.F)


def solve(d, key, x0, x1, N_goal, log_x):
    """secant on ln N against x (or ln x) to the goal row count; returns the iterates"""
    tr = []
    xs = [x0, x1]
    ys = []
    for x in xs:
        N, F = plan(d, **{key: x})
        tr.append(dict(x=x, N=N, F=F))
        ys.append(math.log(N))
        print(key, tr[-1], flush=True)
    for _ in range(6):
        t = (lambda v: math.log(v)) if log_x else (lambda v: v)
        a, b = t(xs[-2]), t(xs[-1])
        nxt = b + (math.log(N_goal) - ys[-1]) * (b - a) / (ys[-1] - ys[-2])
        x = math.exp(nxt) if log_x else nxt
        N, F = plan(d, **{key: x})
        tr.append(dict(x=x, N=N, F=F))
        print(key, tr[-1], flush=True)
        xs.append(x)
        ys.append(math.log(N))
        if abs(N - N_goal) <= 2:
            break
    return tr


for name in NAMES:
    rec = json.load(open(os.path.join(RUNS, name, name + '.json')))
    dsg, val = rec['design'], rec['validation']
    turns = dsg['turns']
    per_leg = [turns[0]] + [turns[i] - turns[i - 1] for i in range(1, len(turns))] + [dsg['N'] - turns[-1]]
    legs_mm = val['length_budget']['leg_path_lengths_mm']
    full = per_leg[1:-1]
    ratio = fig['rows_per_full_leg'] / (sum(full) / len(full))
    d = SP.load(os.path.join(R, 'presets', name + '.yaml'))
    F0, D0 = d['targets']['filtrate_fraction'], d['targets']['cutoff_um']
    N0, Fp = plan(d)
    goal = N0 * ratio
    row = dict(record_code=rec.get('code'), N=dsg['N'], NM=dsg['NM'], rows_per_leg=per_leg,
               full_leg_path_mm=legs_mm[1:-1], footprint_mm=val['drawn_parameters']['footprint_mm'],
               full_leg_rows_vs_fig3a=sum(full) / len(full) / fig['rows_per_full_leg'],
               footprint_length_vs_fig3a_mm=val['drawn_parameters']['footprint_mm'][0] - fig['extent_mm'][0],
               fig3a_over_record_rows=ratio, plan_N_at_targets=N0, plan_F=Fp, F_target=F0, D_um=D0,
               plan_N_goal=goal)
    print(name, json.dumps(row), flush=True)
    tF = solve(d, 'F', F0, F0 + 0.005, goal, False)
    tD = solve(d, 'D', D0, D0 * 0.97, goal, True)
    row['F_to_fill_fig3a_at_D'] = dict(F=tF[-1]['x'], N=tF[-1]['N'], iterates=tF,
                                       filtrate_to_retentate=tF[-1]['x'] / (1.0 - tF[-1]['x']))
    row['D_to_fill_fig3a_at_F'] = dict(D_um=tD[-1]['x'], N=tD[-1]['N'], iterates=tD, D_over_preset=tD[-1]['x'] / D0)
    out['presets'][name] = row
    json.dump(out, open(OUT, 'w'), indent=1)
json.dump(out, open(OUT, 'w'), indent=1)
print('wrote', OUT)
