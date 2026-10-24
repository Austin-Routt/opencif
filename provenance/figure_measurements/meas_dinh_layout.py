"""meas_dinh_layout.py -- Dinh 2024 Fig. 2A main panel ('fig2 large top'; pdfimages dinh-001.png, 982 x 1940): the
whole chip as drawn.  The lateral scale is the chip's width: the arms span the 3.2 mm chip (p. 916), as the inset
shows (arm lanes at the chip's edge).  The inlet port is drawn as a circle, so its height/width tests whether the
panel is scaled the same along and across the chip.  The meander is drawn hatched (grey on the axis), the linear
section solid black.  Writes dinh_layout_meas.json."""
import json
import numpy as np
from PIL import Image

IM = np.asarray(Image.open('../r3/papers/dinh-001.png').convert('L')).astype(float)
AXC = (382, 389)                       # the stem's columns
dark = IM < 128

# inlet port: the dark blob round the axis above the neck (rows 40-112, columns 340-432)
# (columns 350-425 only: label lines reach into the panel at columns 340-348)
rows = [y for y in range(40, 112) if dark[y, 350:425].sum() > 12]
cols = np.nonzero(dark[min(rows):max(rows) + 1, 350:425].any(axis=0))[0] + 350
port_h = max(rows) - min(rows) + 1
port_w = cols.max() - cols.min() + 1
port_cy = 0.5 * (min(rows) + max(rows))
# the arms' span (the chip's width) and the head's top: the widest dark extent in the arm rows
ext = {}
for y in range(180, 300):
    d = np.nonzero(dark[y, 330:440])[0]
    if len(d):
        ext[y] = (d.min() + 330, d.max() + 330)
span = max(b - a + 1 for a, b in ext.values())
head_top = min(y for y, (a, b) in ext.items() if b - a + 1 >= 0.9 * span)
# meander end (hatched -> solid on the axis), outlet start (solid -> open), the drawing's last row
m = np.convolve(IM[:, AXC[0]:AXC[1]].mean(axis=1), np.ones(9) / 9.0, mode='same')
mend = next(y for y in range(900, 1300) if m[y] < 30 and m[y + 20] < 30)
ostart = next(y for y in range(mend + 100, 1900) if m[y] > 60)
last = max(y for y in range(1700, IM.shape[0]) if dark[y, 330:440].any())
SC = 3200.0 / span
# outlet (rows 1800-1846): the middle-channel port is the dark disc on the axis; the two side leads wrap round it.
# Widths by integrated darkness (anti-aliased edges count by their grey level), over a window round each dark run.
def runs(b):
    r, i = [], 0
    while i < len(b):
        if b[i]:
            j = i
            while j < len(b) and b[j]:
                j += 1
            r.append((i, j))
            i = j
        else:
            i += 1
    return r
lead_w, cport_w = [], []
for y in range(1810, 1835):
    row = IM[y]
    rr_ = runs(row < 128)
    blob = next((a, b) for a, b in rr_ if a <= 385 < b)
    cport_w.append(sum((255.0 - row[x]) / 255.0 for x in range(blob[0] - 1, blob[1] + 1)))
    left = max((r for r in rr_ if r[1] <= blob[0] - 2), key=lambda r: r[1])
    right = min((r for r in rr_ if r[0] >= blob[1] + 2), key=lambda r: r[0])
    for a, b in (left, right):
        lead_w.append(sum((255.0 - row[x]) / 255.0 for x in range(a - 1, b + 1)))
port_rows = [y for y in range(1790, 1850) if (IM[y, 380:390] < 128).all() and (IM[y, 365:405] < 128).all()]
# the two ports' centres on the axis (dark runs of the axis column), and the teardrop loop's outer half-width
axcol = IM[1785:1930, 383:388].mean(axis=1) < 128
ax_runs = [(a + 1785, b + 1785) for a, b in runs(axcol) if b - a >= 15]
port_c_row = 0.5 * (ax_runs[0][0] + ax_runs[0][1] - 1)
port_f_row = 0.5 * (ax_runs[-1][0] + ax_runs[-1][1] - 1)
half_w = []
for y in range(1812, 1833):                      # the loop's widest part; columns clear of the label lines
    d = np.nonzero(IM[y, 338:434] < 128)[0]
    if len(d):
        half_w.append(0.5 * (d.max() - d.min() + 1))
out = dict(lateral_scale_um_per_px=SC, arm_span_px=int(span),
           port_px=dict(width=int(port_w), height=int(port_h), centre_row=port_cy),
           port_aspect_h_over_w=port_h / port_w, port_r_um=0.5 * port_w * SC,
           head_top_row=int(head_top), meander_end_row=int(mend), outlet_start_row=int(ostart), last_row=int(last),
           port_centre_to_head_top_mm=(head_top - port_cy) * SC / 1000.0,
           port_centre_to_row1_mm=(head_top - port_cy) * SC / 1000.0
           + json.load(open('dinh_inset_meas.json'))['arm_top_outer_edge_mm_above_row1'][0],
           meander_mm=(mend - head_top) * SC / 1000.0, linear_mm=(ostart - mend) * SC / 1000.0,
           drawing_length_mm=(last - min(rows)) * SC / 1000.0,
           meander_fraction_of_head_to_outlet=(mend - head_top) / float(ostart - head_top),
           meander_fraction_of_drawing=(mend - head_top) / float(last - min(rows)),
           outlet=dict(side_lead_width_px=float(np.median(lead_w)), side_lead_width_um=float(np.median(lead_w)) * SC,
                       side_lead_width_range_um=[float(min(lead_w)) * SC, float(max(lead_w)) * SC],
                       centre_port_width_px=float(max(cport_w)), centre_port_r_um=0.5 * float(max(cport_w)) * SC,
                       centre_port_rows=[int(min(port_rows)), int(max(port_rows))] if port_rows else None,
                       side_lead_width_iqr_um=[float(np.percentile(lead_w, 25)) * SC, float(np.percentile(lead_w, 75)) * SC],
                       port_centre_spacing_um=(port_f_row - port_c_row) * SC,
                       loop_outer_half_width_um=float(np.median(half_w)) * SC,
                       note='rows 1810-1834; widths by integrated darkness; the lead is drawn ~6-7 px at 32 um/px; a few rows '
                            'pick a fragment of the crown of lines entering the middle port (the low end of the range)'),
           note='lateral scale from the arms spanning the 3.2 mm chip; lengths at that scale (the round port shows '
                'the panel is not stretched or squeezed along the chip)')
json.dump(out, open('dinh_layout_meas.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
