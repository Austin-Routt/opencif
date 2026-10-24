"""rc5 (after the rc5 reviews): regenerate every entry of opencif's stretch-end tables from scratch with expfem.py and
compare it with the values shipped in cifsched (STEP_ONE_DEEP, STEP_SYM_DEEP, FLARE_DEEP, FLARE_SHALLOW) and with the
conformal closed form (the shallow steps).  The meshes are the ones the tables were made with: steps hmin 0.005,
hmax 0.08 x max(1, W/3); flares hmin 0.005, hmax 0.08; the one-sided R 2.4 entry on the finer mesh (hmin 0.0025,
hmax 0.06) its table value came from.  One worker by default (WORKERS=2 for two: a first attempt with two workers
beside another solve ran out of the 6 GB memory limit, the larger symmetric steps taking 3.5 GB each); each result is
printed as it finishes, in table order.  About 15 min on one core (53 entries, 14.3 min here).
Usage: [WORKERS=1] python3 regen_fem_tables.py OPENCIF_ROOT [OUT.json]   (needs requirements-provenance.txt)"""
import json, math, os, sys, time
from multiprocessing import Pool
R = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else 'regen_fem_tables.json'
sys.path[:0] = [os.path.dirname(os.path.abspath(__file__)), os.path.join(R, 'opencif', 'engine')]
import cifsched as S


def job(spec):
    import expfem as E
    kind, W, deg = spec
    t = time.time()
    if kind in ('one', 'sym'):
        if kind == 'one' and abs(W - 2.4) < 1e-9:
            hmin, hmax = 0.0025, 0.06
        else:
            hmin, hmax = 0.005, 0.08 * max(1.0, W / 3.0)
        d, ne = E.excess('step', W, hmin=hmin, hmax=hmax, limit='deep', sym=(kind == 'sym'))
        s, _ = E.excess('step', W, hmin=hmin, hmax=hmax, limit='shallow', sym=(kind == 'sym'))
    else:
        d, ne = E.excess('flare', W, deg, hmin=0.005, hmax=0.08, limit='deep')
        s, _ = E.excess('flare', W, deg, hmin=0.005, hmax=0.08, limit='shallow')
    return dict(kind=kind, ratio=W, deg=deg, deep=d, shallow=s, elems=ne, seconds=round(time.time() - t, 1))


if __name__ == '__main__':
    jobs = [('one', W, None) for W, _ in S.STEP_ONE_DEEP if W > 1] + [('sym', W, None) for W, _ in S.STEP_SYM_DEEP if W > 1]
    jobs += [('flare', W, g) for W in S.FLARE_RATIOS if W > 1 for g in S.FLARE_DEGS if g > 0]
    res = []
    with Pool(int(os.environ.get('WORKERS', '1')), maxtasksperchild=1) as p:
        for r in p.imap(job, jobs):
            res.append(r)
            print('%-5s R %.2f deg %-5s deep %.5f shallow %.5f elems %6d (%.0f s)' % (r['kind'], r['ratio'], r['deg'], r['deep'],
                  r['shallow'], r['elems'], r['seconds']), flush=True)
    json.dump(res, open(OUT, 'w'), indent=1)
    worst = {}
    for r in res:
        if r['kind'] == 'flare':
            i, j = S.FLARE_RATIOS.index(r['ratio']), S.FLARE_DEGS.index(r['deg'])
            tab = [(S.FLARE_DEEP[i][j], r['deep']), (S.FLARE_SHALLOW[i][j], r['shallow'])]
            key = 'flare (deep and shallow tables)'
        else:
            table = dict(S.STEP_ONE_DEEP if r['kind'] == 'one' else S.STEP_SYM_DEEP)
            tab = [(table[r['ratio']], r['deep'])]
            key = '%s-step deep table' % r['kind']
            cf = S.step_shallow(r['ratio'], r['kind'] == 'one')
            worst['%s-step shallow: closed form vs regenerated FEM' % r['kind']] = max(
                worst.get('%s-step shallow: closed form vs regenerated FEM' % r['kind'], 0.0), abs(cf - r['shallow']))
        worst[key] = max(worst.get(key, 0.0), max(abs(a - b) for a, b in tab))
    for k, v in worst.items():
        print('largest |shipped - regenerated| %-50s %.1e' % (k, v))
