"""test_lattice.py: the option lattice at the plan stage (review FINAL_REVIEW, section 5.3).

The lattice crosses every named process with every outlet style and with F at both extremes and the middle, on a
paperclip fitted to a 25 x 75 mm chip, and adds the other layouts at F 0.9 and a center channel at the plateau edge.
Each case is planned and gated without drawing (generate(check_only='plan')), which includes the outlet's own plan on
the leg counts that the fitter would try (pipeline.outlet_precheck). The full check_only=True also draws the outline in
memory and checks its collisions, which tests/test_composition.py covers.

For every case, the answer must be a record with exit code 0, 3, or 5, or a SpecError (exit 6), and nothing else may
be raised. No message may carry a library's text (a root finder's "f(a) and f(b) must have different signs", or a
traceback), and an outlet refusal must name the outlet. The table is printed. Run it with
`python3 tests/test_lattice.py` (10 to 20 minutes for about 40 plans, each with its outlet plan on up to ten leg
counts)."""
import os, sys, time, traceback
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import tempfile as _tempfile
# Scratch output goes outside the package tree, so that a read-only install works, and OPENCIF_TEST_OUT overrides
# the location (review R3 of rc4).
T_OUT = os.environ.get('OPENCIF_TEST_OUT') or _tempfile.mkdtemp(prefix='opencif_test_')
if not os.environ.get('OPENCIF_TEST_OUT'):               # a temporary directory of our own, removed at exit
    import atexit as _atexit, shutil as _shutil
    _atexit.register(_shutil.rmtree, T_OUT, True)

sys.path.insert(0, ROOT)
import opencif
from opencif import spec as SP
from opencif import pipeline as _PL_CODE                # noqa: E402  (every suite's log names the code it ran on)
print('opencif %s (code %s)' % (_PL_CODE.CODE['version'], _PL_CODE.CODE['sha256']))

OUT = T_OUT
LIBRARY_TEXT = ('f(a) and f(b)', 'different signs', 'Traceback', 'brentq', 'scipy', 'numpy')
D_OF = {'generic19': 4.0, 'gifford2014': 5.2, 'dinh2024': 2.5, 'xia2016': 7.0}
CASES = []
for proc, D in D_OF.items():
    for outlet in ('separate', 'teardrop', 'ring'):
        for F in (0.1, 0.5, 0.9):
            CASES.append(('%s %s F%.1f paperclip' % (proc, outlet, F),
                          dict(cutoff_um=D, filtrate_fraction=F, process=proc, chip=(25, 75), layout='paperclip',
                               outlet=outlet, mode='evaluate')))
for layout in ('linear', 'serpentine'):
    for outlet in ('separate', 'teardrop'):
        CASES.append(('generic19 %s F0.9 %s' % (outlet, layout),
                      dict(cutoff_um=4.0, filtrate_fraction=0.9, chip=(25, 75), layout=layout, outlet=outlet)))
CASES.append(('generic19 separate F0.5 linear, u = 0.035 (plateau edge)',
              dict(cutoff_um=4.0, filtrate_fraction=0.5, width_in_um=4.0 / 0.035, layout='linear')))
CASES.append(('generic19 ring F0.9 linear', dict(cutoff_um=4.0, filtrate_fraction=0.9, layout='linear', outlet='ring')))

fails, rows = [], []
for name, kw in CASES:
    t0 = time.time()
    try:
        r = opencif.generate(**kw, check_only='plan', out=OUT)
        code, msg = r.exit_code, r.verdict
        ok = code in (0, 3, 5) and not any(t in msg for t in LIBRARY_TEXT)
        if code == 5 and (r.record.get('outlet_refusal') or ''):
            ok = ok and kw.get('outlet', 'separate') in r.record['outlet_refusal']
    except SP.SpecError as e:
        code, msg = 6, str(e)
        ok = not any(t in msg for t in LIBRARY_TEXT)
    except Exception as e:                                       # anything else is a defect
        code, msg, ok = 'EXC', '%s: %s | %s' % (type(e).__name__, e, traceback.format_exc()[-300:]), False
    rows.append((name, code, round(time.time() - t0, 1), msg))
    print('%-58s %-4s %5.1f s  %s %s' % (name, code, time.time() - t0, 'ok ' if ok else 'FAIL', msg[:150]), flush=True)
    if not ok:
        fails.append(name)
from collections import Counter
print('%d cases; exit codes %s' % (len(rows), dict(Counter(str(r_[1]) for r_ in rows))))
print('ALL OK' if not fails else 'FAILED: %s' % fails)
sys.exit(1 if fails else 0)
