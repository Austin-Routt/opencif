"""Run the opencif CLI from a scratch tree with outlet_exp_domain's table end forced to a lower ratio (test only).
Usage: python3 run_forced.py TREE R_END|none ARGS..."""
import os, runpy, sys
tree, r_end, args = sys.argv[1], sys.argv[2], sys.argv[3:]
sys.path.insert(0, tree)
os.chdir(tree)
import opencif  # noqa: F401  (puts the engine on the path)
import cifsched as S
assert os.path.abspath(S.__file__).startswith(os.path.abspath(tree)), S.__file__
calls = []
if r_end != 'none':
    R_END = float(r_end)

    def forced(spec, w_a, w_b):
        wn, ww = min(w_a, w_b), max(w_a, w_b)
        if getattr(spec, 'junction_band', None) == 'none' or ww <= wn * (1.0 + 1e-6):
            return None
        R = ww / wn
        calls.append(R)
        return ('forced test: width ratio %.3f beyond %.2f' % (R, R_END)) if R > R_END else None
    S.outlet_exp_domain = forced
sys.argv = ['opencif'] + args
try:
    runpy.run_module('opencif', run_name='__main__', alter_sys=True)
finally:
    if calls:
        print('[test] outlet_exp_domain called %d times; ratios %s' % (len(calls), sorted(set(round(c, 3) for c in calls))[:8]))
