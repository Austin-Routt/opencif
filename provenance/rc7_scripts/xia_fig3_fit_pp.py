"""rc6: Xia RS-7 with BOTH decoupled targets read off Fig. 3 instead of the text (the rc6 traces report, sec. 12.4).

F from panel g: the last side channel's width on Fig. 3g (results/xia/xia_results_rc6.json, measured.left_side_channel_um,
outline to outline at the rc6 scale) is what the schedule's last straight side channel must equal; the side width at
the end is set by F (the outlet identity), nearly independently of D.  D from Fig. 3a: the rows per full leg that fill
Fig. 3a's legs at the measured row pitch, as xia_leg_fill.py computes them (its plan_N_goal: the preset's planned N
times Fig. 3a's ratio of full-leg rows to the record's).  The two are solved alternately (secant steps, as
xia_leg_fill.py) until neither moves: F for the width at the current D, D for the row count at that F.

The text's values (printed f_ref -> D 6.8644 um, printed ~15:1 -> F 0.938) stay in the preset; this is a reading of the
figure for the variant variants/xia2016_rs7_depicted_fig3.yaml, and it inherits the model's premises (the frozen
cutoff law's kernel, the plateau prefactor, the 300 -> 150 um constant-shear taper; xia_leg_fill.py states them).

Usage: python3 -I xia_fig3_fit.py OPENCIF_ROOT LEG_FILL.json TRACE_G.json OUT.json"""
import copy
import json
import math
import os
import sys

R, LF, TG, OUT = sys.argv[1:5]
NAME = 'xia2016_rs7_depicted'
sys.path[:0] = [R]
from opencif import spec as SP
from opencif.physics import schedule as S

lf = json.load(open(LF))['presets'][NAME]
N_goal = float(lf['plan_N_goal'])
tg = json.load(open(TG))
w_goal = float(tg['measured']['left_side_channel_um'])
d0 = SP.load(os.path.join(R, 'presets', NAME + '.yaml'))
d0.setdefault('centre_channel', {})['shear_formula'] = 'parallel_plates'   # rc7 (XIA-1)
F0, D0 = d0['targets']['filtrate_fraction'], d0['targets']['cutoff_um']
log = []


def plan(F, D):
    d = copy.deepcopy(d0)
    d['targets'].update(filtrate_fraction=F, cutoff_um=D)
    sp, _ = SP.to_engine(SP.load(d))
    sc = S.plan(sp)
    row = dict(F=F, D_um=D, N=int(sc.N), ws_last_um=float(sc.ws[sc.N - 2]), F_planned=float(sc.F))
    log.append(row)
    print(json.dumps(row), flush=True)
    return row


def secant(f, x0, x1, goal, tol, n=6):
    a, b = x0, x1
    ya, yb = f(a), f(b)
    for _ in range(n):
        if abs(yb - goal) <= tol or yb == ya:
            break
        a, b, ya = b, b + (goal - yb) * (b - a) / (yb - ya), yb
        yb = f(b)
    return b, yb


F, D = F0, float(lf['D_to_fill_fig3a_at_F']['D_um']) if 'D_to_fill_fig3a_at_F' in lf else D0
for it in range(3):
    F_new, w = secant(lambda x: plan(x, D)['ws_last_um'], F, F + 0.002, w_goal, 0.3)
    D_new, lnN = secant(lambda x: math.log(plan(F_new, math.exp(x))['N']), math.log(D), math.log(D * 0.995),
                        math.log(N_goal), 2.0 / N_goal)
    D_new = math.exp(D_new)
    moved = abs(F_new - F) < 2e-5 and abs(D_new - D) < 2e-4
    F, D = F_new, D_new
    if moved:
        break
fin = plan(F, D)
out = dict(method=__doc__.split('\n\n')[1].replace('\n', ' '), preset=NAME,
           goals=dict(last_side_channel_um=w_goal, from_trace=os.path.basename(TG), N=N_goal, from_leg_fill=os.path.basename(LF)),
           text_values=dict(F=F0, D_um=D0, filtrate_to_retentate=F0 / (1 - F0)),
           result=dict(F=F, D_um=D, N=fin['N'], ws_last_um=fin['ws_last_um'], filtrate_to_retentate=F / (1 - F),
                       D_over_text=D / D0),
           iterates=log)
json.dump(out, open(OUT, 'w'), indent=1)
print('result', json.dumps(out['result']))
