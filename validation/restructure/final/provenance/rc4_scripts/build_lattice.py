"""Full builds, on rc4 with two workers, of every plan-stage lattice case that plans inside the envelope (21) and of one
Xia-process case that plans outside it (the ring at F 0.5; the other seven Xia plans are the same M13 exit 3); one JSON
per case and a summary.  Evidence for the report (does every variation that plans also draw and pass?), not a shipped suite.
ROOT and OUT below are the paths of the session that ran it; set them to your tree and an output folder."""
import json, os, sys, time, traceback
from multiprocessing import Pool
ROOT = 'opencif_rc4'
OUT = 'scratch/rc4work/blat'
D_OF = {'generic19': 4.0, 'gifford2014': 5.2, 'dinh2024': 2.5}
CASES = []
for proc, D in D_OF.items():
    for outlet in ('separate', 'teardrop', 'ring'):
        for F in (0.5, 0.9):
            if outlet == 'ring' and F == 0.9 and proc != 'dinh2024':
                continue                                   # refused before drawing (the lattice)
            CASES.append(('%s_%s_F%s_paperclip' % (proc, outlet, str(F).replace('.', 'p')),
                          dict(cutoff_um=D, filtrate_fraction=F, process=proc, chip=(25, 75), layout='paperclip',
                               outlet=outlet, mode='evaluate')))
for layout in ('linear', 'serpentine'):
    for outlet in ('separate', 'teardrop'):
        CASES.append(('generic19_%s_F0p9_%s' % (outlet, layout),
                      dict(cutoff_um=4.0, filtrate_fraction=0.9, chip=(25, 75), layout=layout, outlet=outlet)))
CASES.append(('generic19_separate_F0p5_linear_u0p035', dict(cutoff_um=4.0, filtrate_fraction=0.5, width_in_um=4.0 / 0.035,
                                                            layout='linear')))
CASES.append(('xia2016_ring_F0p5_paperclip', dict(cutoff_um=7.0, filtrate_fraction=0.5, process='xia2016', chip=(25, 75),
                                                  layout='paperclip', outlet='ring', mode='evaluate')))


def run(case):
    name, kw = case
    sys.path.insert(0, ROOT)
    import opencif
    t0 = time.time()
    rec = dict(name=name, kwargs=json.loads(json.dumps(kw, default=list)))
    try:
        r = opencif.generate(**kw, name=name, out=os.path.join(OUT, name), render=False, verbose=False)
        a = (r.record or {}).get('assessment') or {}
        rec.update(exit_code=r.exit_code, verdict=r.verdict, status=a.get('status'),
                   failed=[c['check'] for c in a.get('checks', []) if not c['ok'] and c['kind'] == 'hard'])
    except Exception as e:
        rec.update(exception=type(e).__name__, message=str(e)[:800], trace=traceback.format_exc()[-1500:])
    rec['seconds'] = round(time.time() - t0, 1)
    json.dump(rec, open(os.path.join(OUT, name + '.json'), 'w'), indent=1, default=str)
    return name, rec.get('exit_code', rec.get('exception')), rec['seconds'], (rec.get('verdict') or rec.get('message') or '')[:160]


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    print('%d cases' % len(CASES), flush=True)
    with Pool(2) as p:
        for n, c, dt, v in p.imap_unordered(run, CASES):
            print('%s %-45s %-6s %6.0f s  %s' % (time.strftime('%H:%M:%S'), n, c, dt, v), flush=True)
    print('BLAT_DONE', flush=True)
