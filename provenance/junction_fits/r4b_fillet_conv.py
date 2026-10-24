import sys, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from r4b_fillet import deep, hele_shaw
for rf in (0.0, 0.18):
    for n in (40, 60, 90, 120):
        t0 = time.time()
        print('r_f/w %.2f n %3d: deep c_T %.4f  (%.0fs)' % (rf, n, deep(n, rf, Ld=2.0, Lu=2.0), time.time() - t0), flush=True)
    print('r_f/w %.2f: shallow n160 %.4f' % (rf, hele_shaw(160, rf, Ld=2.0, Lu=2.0)), flush=True)
