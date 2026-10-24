"""Regression of the rc6 Dinh length ports: run dinh_length_bound.py, dinh_taper_alternatives.py and
dinh_schedule_at_F.py on the rc5 post-review tree with the rc5 pitch tolerance (--pitch-tol 1.0; the rc5 presets'
pitch is 19 + 27 = 46 um) and compare with rc5's own outputs (provenance/rc5_scripts/).  The only allowed differences
are the renamed length key (length_mm_at_46 -> length_mm_at_pitch) and the added inputs (preset, pitch_um,
pitch_tolerance_um).  Usage: python3 regress_length_rc5.py RC5_ROOT WORKDIR  (prints the verdicts; exit 1 on a difference)"""
import json
import os
import subprocess
import sys

R5, W = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
H = os.path.dirname(os.path.abspath(__file__))
O5 = os.path.join(R5, 'provenance', 'rc5_scripts')
os.makedirs(W, exist_ok=True)
py = sys.executable
bad = 0


def run(args, out):
    with open(out, 'w') as fh:
        subprocess.run([py] + args, check=True, stdout=fh, cwd=W)


for k in ('depicted', 'specified'):
    j = os.path.join(W, 'bound_%s.json' % k)
    run([os.path.join(H, 'dinh_length_bound.py'), R5, os.path.join(R5, 'presets', 'dinh2024_v1_%s.yaml' % k), j,
         '--pitch-tol', '1.0'], os.path.join(W, 'bound_%s.out' % k))
    a, b = json.load(open(os.path.join(O5, 'dinh_length_bound_%s.json' % k))), json.load(open(j))
    for r in b['bounds']:
        r['length_mm_at_46'] = r.pop('length_mm_at_pitch')
    extra = sorted(x for x in b if x not in a)
    for x in extra:
        b.pop(x)
    ok = a == b
    bad += not ok
    print('dinh_length_bound %-9s %s (added inputs: %s)' % (k, 'identical' if ok else 'DIFFERENT', ', '.join(extra)))
j = os.path.join(W, 'taper.json')
run([os.path.join(H, 'dinh_taper_alternatives.py'), R5, os.path.join(R5, 'presets', 'dinh2024_v1_depicted.yaml'), j,
     os.path.join(W, 'bound_depicted.json'), os.path.join(W, 'bound_specified.json'), '--pitch-tol', '1.0'],
    os.path.join(W, 'taper.out'))
a, b = json.load(open(os.path.join(O5, 'dinh_taper_alternatives.json'))), json.load(open(j))
b.pop('pitch_tolerance_um')
ok = a == b and open(os.path.join(O5, 'dinh_taper_alternatives.out')).read() == open(os.path.join(W, 'taper.out')).read()
bad += not ok
print('dinh_taper_alternatives    %s (JSON and printed output; added input: pitch_tolerance_um)' % ('identical' if ok else 'DIFFERENT'))
j = os.path.join(W, 'schedule.json')
run([os.path.join(H, 'dinh_schedule_at_F.py'), R5, j, '--runs', os.path.join(R5, 'runs')], os.path.join(W, 'schedule.out'))
ok = (json.load(open(os.path.join(O5, 'dinh_schedule_at_F.json'))) == json.load(open(j)) and
      open(os.path.join(O5, 'dinh_schedule_at_F.out')).read() == open(os.path.join(W, 'schedule.out')).read())
bad += not ok
print('dinh_schedule_at_F         %s (JSON and printed output; x_top from the plan = the records\' x_top)'
      % ('identical' if ok else 'DIFFERENT'))
sys.exit(1 if bad else 0)
