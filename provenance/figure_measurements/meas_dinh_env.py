import json
import numpy as np
from PIL import Image
IM = np.asarray(Image.open('../r3/papers/dinh-001.png').convert('L')).astype(float)
SC, XC, Y_ROW1 = 9.89, 192.6, 599.5
env = []
for y in range(185, 599):
    row = IM[y, 110:276]
    blk = np.nonzero(row < 128)[0]
    if len(blk) == 0:
        continue
    a, b = blk[0], blk[-1]
    def edge(i, j):
        vi, vj = row[i], row[j]
        return i + (128 - vi) / (vj - vi) * (j - i)
    xl = edge(a, a - 1) + 110
    xr = edge(b, b + 1) + 110
    env.append(((y - Y_ROW1) * SC, (xr - xl) / 2 * SC, ((xl + xr) / 2 - XC) * SC))
env = np.array(env)
d = json.load(open('dinh_filter_meas.json'))
d['envelope'] = env.tolist()
json.dump(d, open('dinh_filter_meas.json', 'w'), indent=1)
for r in env[::8]:
    print('x %8.1f half %6.1f c %5.1f' % tuple(r))
