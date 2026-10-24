import json
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, gaussian_filter1d
IM = np.asarray(Image.open('../r3/papers/xia-002.png').convert('L')).astype(float)
HP = IM - gaussian_filter(IM, 6.0)
SC = 250.0 / 128.0
XC = 242.0
res = []
for k in range(11):
    y0 = 100 + int(round(65.3 * k)); y1 = min(y0 + 66, 752)
    if y1 - y0 < 30: break
    cnt = (HP[y0:y1, 400:540] < -8).sum(1)
    # horizontal outline rows: local maxima of the count
    ys = [y0 + i for i in range(1, len(cnt) - 1) if cnt[i] >= cnt[i - 1] and cnt[i] > cnt[i + 1] and cnt[i] > 80]
    ends = []
    for y in ys:
        line = gaussian_filter1d(HP[y - 1:y + 2, :].min(0), 0.7)
        # walk right from x=420 while the line stays dark (allow gaps <= 4 px)
        x, gap = 420, 0
        while x < 735 and gap <= 4:
            gap = 0 if line[x] < -6 else gap + 1
            x += 1
        ends.append((y, x - gap - 1))
    res.append(dict(row=k + 1, lines=ends))
    print(k + 1, ends)
json.dump(res, open('xia_head_lines.json', 'w'), indent=1)
