"""Print the ring outlet's side-lead resistance target while a build plans its layout (why rc3 failed the ring at F 0.9).
usage: diag_ring2.py ROOT '{generate kwargs as JSON}' OUTDIR"""
import sys, json
sys.path.insert(0, sys.argv[1])
import opencif, opencif.pipeline
L = sys.modules['layout']; S = sys.modules['cifsched']
orig_w = L.ring2_width; seen = set()
def spy(g, target, h, spec=None):
    k = round(target, 3)
    if k not in seen:
        seen.add(k)
        print('ring2_width: target %.4e  R(w=1e5) %.4e  R(w=40) %.4e  L %.1f um  rq %.1f' % (target, S.rr(1.0e5, h) * g['L'] * S.UM, S.rr(40.0, h) * g['L'] * S.UM, g['L'], g['rq']), file=sys.stderr)
    return orig_w(g, target, h, spec)
L.ring2_width = spy
orig_cp = L.combined_plan
def cp(sc):
    spec, h, N = sc.spec, sc.spec.h, sc.N
    w_last = sc.ws[N - 2]; r_str = S.rr(w_last, h) * spec.L_out * S.UM
    print('combined_plan: N %d w_last %.2f wc_end %.1f Rs_out %.4e slope*? P/S %.4e r_str %.4e R_ext_side %.4e L_out %.1f F %.3f' % (
        N, w_last, sc.wc[-1], sc.Rs_out, sc.P[N-1]/sc.S[N-1], r_str, L.R_ext_side(spec), spec.L_out, getattr(spec,'F',float('nan'))), file=sys.stderr)
    return orig_cp(sc)
L.combined_plan = cp
kw = json.loads(sys.argv[2]); kw['chip'] = tuple(kw['chip'])
try:
    r = opencif.generate(**kw, check_only=False, out=sys.argv[3], verbose=False, render=False)
    print(r.exit_code, r.verdict[:300])
except Exception as e:
    print('EXC', type(e).__name__, str(e)[:200])
