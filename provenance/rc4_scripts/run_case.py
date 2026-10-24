"""Run one generator case on a given opencif tree. usage: run_case.py ROOT OUTDIR NAME CHECK_ONLY JSON_KWARGS"""
import sys, os, json, time, traceback
root, outdir, name, co, kw = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] == '1', json.loads(sys.argv[5])
sys.path.insert(0, root)
import opencif
for k in ('chip',):
    if k in kw and isinstance(kw[k], list): kw[k] = tuple(kw[k])
if 'taper' in kw and isinstance(kw['taper'], list): kw['taper'] = tuple(kw['taper'])
t0 = time.time(); rec = dict(name=name, kwargs=kw, check_only=co)
try:
    r = opencif.generate(**kw, check_only=co, out=os.path.join(outdir, name), verbose=False, render=False)
    rec.update(exit_code=r.exit_code, verdict=r.verdict, ok=r.ok)
    a = (r.record or {}).get('assessment') or {}
    rec['status'] = a.get('status'); rec['notes'] = a.get('notes') or (r.record or {}).get('notes')
except Exception as e:
    rec.update(exception=type(e).__name__, message=str(e)[:1500], trace=traceback.format_exc()[-3000:])
rec['seconds'] = round(time.time() - t0, 1)
json.dump(rec, open(os.path.join(outdir, name + '.json'), 'w'), indent=1, default=str)
print(name, rec.get('exit_code', rec.get('exception')), rec['seconds'], 's'); print((rec.get('verdict') or rec.get('message') or '')[:800])
