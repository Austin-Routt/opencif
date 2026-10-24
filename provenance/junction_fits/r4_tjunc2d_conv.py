import sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from r4_tjunc2d import deep, hele_shaw
for stem in (0.0, 2.0):
    v = [deep(n, Ld=3.0, Lu=3.0, stem=stem) for n in (20, 40, 80)]
    p = __import__('math').log(abs((v[0] - v[1]) / (v[1] - v[2]))) / __import__('math').log(2)
    print('deep stem %.0f: n 20/40/80 %s, observed order %.2f, extrapolated %.4f' % (stem, ['%.4f' % x for x in v], p, v[2] + (v[2] - v[1]) / (2 ** p - 1)), flush=True)
    s = [hele_shaw(n, Ld=3.0, Lu=3.0, stem=stem) for n in (20, 40, 80)]
    print('shallow stem %.0f: %s' % (stem, ['%.4f' % x for x in s]), flush=True)
