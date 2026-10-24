"""rc5 (after the second independent check): the flared outlet's stretch rules on layouts no shipped spec builds.
Adds outlet.flare_deg 14.5 to three examples (a linear separate outlet whose lead does not widen; a teardrop on a
serpentine and on a paperclip, at one and three legs), plans each layout and prints the stretch the layout sets:
no widening keeps the spec's stretch; one leg (no U-turn before the outlet) gives flare_start + the flare + its end arc;
three legs give the loop's clearance stretch + the flare + its end arc.  Plans only (no drawing).
Usage: python3 flare_rules_check.py OPENCIF_ROOT"""
import copy, os, sys
R = sys.argv[1]
sys.path[:0] = [R, os.path.join(R, 'opencif'), os.path.join(R, 'opencif', 'engine')]
from opencif import spec as SP
import cifsched as S, layout as L
CASES = [('m01_linear_const_pack_plain_separate', 'linear', 1), ('m05_serpentine_const_arms_plain_teardrop', 'serpentine', 3),
         ('m05_serpentine_const_arms_plain_teardrop', 'serpentine', 1), ('m09_paperclip_const_narrow_plain_teardrop', 'paperclip', 3),
         ('m09_paperclip_const_narrow_plain_teardrop', 'paperclip', 1)]
for name, kind, legs in CASES:
    d = copy.deepcopy(SP.load(os.path.join(R, 'examples', 'matrix', name + '.yaml')))
    d.setdefault('outlet', {})['flare_deg'] = 14.5
    sp, op = SP.to_engine(SP.load(d))
    sc = S.plan(sp)
    sc2, _ = L.plan_layout(sc, kind, legs, None)
    N, i = sc2.N, sc2.info
    Lf = i.get('flare_len', 0.0)
    tail = Lf + (L.flare_kink(sc2.spec, Lf)[1] if Lf > 0 else 0.0)
    print('%-44s %-10s legs %d  w_s(N) %6.2f  lead %7.2f  rule %-9s  base %8.2f  flare %6.2f  stretch %8.2f  (base or 67, + flare + arc: %8.2f)'
          % (name, kind, legs, sc2.ws[N - 2], i.get('w_lead'), i.get('flare_kind'), i.get('L_out_base'), Lf, sc2.spec.L_out,
             (i.get('L_out_base') if Lf <= 0 else (max(i.get('L_out_base'), 67.0) if i.get('flare_kind') == 'clearance' else 67.0) + tail)))
