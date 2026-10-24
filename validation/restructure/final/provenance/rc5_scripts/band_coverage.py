"""rc5 (after the rc5 reviews): what the stretch-end term's model band (cifsched: the depth weight moved by
+-max(weight, 0.1) and c scaled by -+30 %, junction_band 'lo' / 'hi') spans on every shipped design that has the term,
against two things it stands in for:
  * the depth proxy itself: c lies between its deep and shallow limits whatever the weight; is that whole range inside
    the band?
  * the teardrop lead's first quarter arc (centreline radius 2 lead widths), which starts at the stretch's end: the 2-D
    interaction of step_bend.py (coarse mesh, interpolated in ln R between the solved ratios 2.4, 2.9 and 3.87) added
    to both limits of a one-sided step; is the shifted c inside the band?  (Ring leads turn on their own radius and
    angle and the flares end in a different corner: not computed, reported as such.)
Reads every record's design.outlet.stretch_end.  After an independent check of the rc5-review report: the first
question is also asked in general, over width ratio and depth weight, for the one-sided and centred steps (cifsched's
own tables and closed form), since containing both limits on the shipped designs is not a property of the band.
Usage: python3 band_coverage.py OPENCIF_ROOT [step_bend.out]"""
import glob, json, math, os, re, sys
R = sys.argv[1]
sys.path[:0] = [os.path.join(R, 'opencif', 'engine')]
import cifsched as CS
SB = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'step_bend.out')
WIDEN = 0.3                                   # cifsched.OUTLET_EXP_WIDEN

inter = {}
for line in open(SB):
    m = re.search(r'mesh hmin ([\d.]+) hmax ([\d.]+)\s+W ([\d.]+).*deep ([-+][\d.]+) \(([-+][\d.]+) % of c\) shallow '
                  r'([-+][\d.]+) \(([-+][\d.]+) % of c\)', line)
    if m and abs(float(m.group(1)) - 0.005) < 1e-9:
        inter[float(m.group(3))] = (float(m.group(5)) / 100.0, float(m.group(7)) / 100.0)
Ws = sorted(inter)


def interaction(Rr):
    if not Ws or Rr < Ws[0] - 0.01 or Rr > Ws[-1] + 0.01:   # solved ratios are rounded (wbc's 3.8702 is solved as 3.87)
        return None
    Rr = min(max(Rr, Ws[0]), Ws[-1])
    lr = math.log(Rr)
    for a, b in zip(Ws, Ws[1:] + [Ws[-1]]):
        if a - 1e-9 <= Rr <= b + 1e-9:
            t = 0.0 if b == a else (lr - math.log(a)) / (math.log(b) - math.log(a))
            return tuple(inter[a][k] + t * (inter[b][k] - inter[a][k]) for k in (0, 1))


rows = []
for jp in sorted(glob.glob(os.path.join(R, 'runs', '*', '*.json'))):
    nm = os.path.basename(os.path.dirname(jp))
    if os.path.basename(jp) != nm + '.json':
        continue
    o = (json.load(open(jp)).get('design') or {}).get('outlet', {})
    se = o.get('stretch_end')
    if not se or se.get('kind') == 'none':
        continue
    lam, dp, sh, c = se['depth_weight'], se['c_deep'], se['c_shallow'], se['c']
    l_lo, l_hi = max(lam - max(lam, 0.1), 0.0), min(lam + max(lam, 0.1), 1.0)
    c_lo = (dp + l_lo * (sh - dp)) * (1.0 + WIDEN)
    c_hi = (dp + l_hi * (sh - dp)) * (1.0 - WIDEN)
    b0, b1 = min(c_lo, c_hi), max(c_lo, c_hi)
    row = dict(design=nm, style=o.get('style'), kind=se['kind'], ratio=round(se['ratio'], 3), weight=round(lam, 3),
               c=round(c, 4), deep_rel=round(dp / c - 1.0, 3), shallow_rel=round(sh / c - 1.0, 3),
               band_rel=[round(b0 / c - 1.0, 3), round(b1 / c - 1.0, 3)],
               limits_inside_band=bool(b0 <= min(dp, sh) + 1e-12 and max(dp, sh) <= b1 + 1e-12))
    it = interaction(se['ratio']) if (se['kind'] == 'step, one-sided' and o.get('style') == 'teardrop') else None
    if it is not None:
        cb = dp * (1.0 + it[0]) + lam * (sh * (1.0 + it[1]) - dp * (1.0 + it[0]))
        row.update(arc_interaction_deep_shallow=[round(it[0], 3), round(it[1], 3)], c_with_arc_rel=round(cb / c - 1.0, 3),
                   c_with_arc_inside_band=bool(b0 <= cb <= b1))
    else:
        row.update(arc_interaction='not computed (%s)' % ('ring lead: its own turn radius and angle' if o.get('style') == 'ring'
                                                           else se['kind']))
    rows.append(row)
    print(json.dumps(row), flush=True)
summ = dict(designs=len(rows), interaction_solved_at=Ws,
            limits_inside_band=sum(r['limits_inside_band'] for r in rows),
            teardrop_steps=[r['design'] for r in rows if 'c_with_arc_rel' in r],
            c_with_arc_rel_range=[min(r['c_with_arc_rel'] for r in rows if 'c_with_arc_rel' in r),
                                  max(r['c_with_arc_rel'] for r in rows if 'c_with_arc_rel' in r)],
            c_with_arc_inside_band=sum(r.get('c_with_arc_inside_band', False) for r in rows))
print('SUMMARY', json.dumps(summ))


def band_of(dp, sh, lam):
    l_lo, l_hi = max(lam - max(lam, 0.1), 0.0), min(lam + max(lam, 0.1), 1.0)
    a, b = (dp + l_lo * (sh - dp)) * (1.0 + WIDEN), (dp + l_hi * (sh - dp)) * (1.0 - WIDEN)
    return min(a, b), max(a, b)


general = {}
for kind, tab, one in (('step, one-sided', CS.STEP_ONE_DEEP, True), ('step, centred', CS.STEP_SYM_DEEP, False)):
    out = []
    for lam100 in range(0, 101, 2):
        lam = lam100 / 100.0
        runs, cur = [], None
        for k in range(0, 201):
            Rr = math.exp(math.log(1.05) + k * (math.log(8.0) - math.log(1.05)) / 200.0)
            dp, sh = CS._interp_logR(tab, Rr), CS.step_shallow(Rr, one)
            b0, b1 = band_of(dp, sh, lam)
            out_ = not (b0 <= min(dp, sh) + 1e-12 and max(dp, sh) <= b1 + 1e-12)
            if out_ and cur is None:
                cur = [Rr, Rr]
            elif out_:
                cur[1] = Rr
            elif cur is not None:
                runs.append(cur); cur = None
        if cur is not None:
            runs.append(cur)
        out.append(dict(weight=lam, ratios_where_a_limit_is_outside_the_band=[[round(a_, 2), round(b_, 2)] for a_, b_ in runs]))
    general[kind] = out
    print('GENERAL %s (ratio 1.05-8, weight 0-1 by 0.02): a limit outside the band at' % kind)
    for o in out:
        if o['ratios_where_a_limit_is_outside_the_band']:
            print('   weight %.2f: R %s' % (o['weight'], ', '.join('%.2f-%.2f' % tuple(r) for r in o['ratios_where_a_limit_is_outside_the_band'])))
json.dump(dict(rows=rows, summary=summ, general=general), open('band_coverage.json', 'w'), indent=1)
