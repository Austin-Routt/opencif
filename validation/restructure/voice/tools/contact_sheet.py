"""contact_sheet.py: a one-look review of every build in a single PDF.

Each record's layout figure is grouped as in docs/RESULTS.md and titled with its verdict, the drawn components (head
style and arms, sealing, filter posts, outlet, and legs), the no-intersection check, and the check that the drawn
components are those requested.

Usage: python3 tools/contact_sheet.py [OUT.pdf]   (default docs/CONTACT_SHEET.pdf)"""
import glob
import json
import os
import sys
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GROUPS = [('Reproductions (presets/)', 'presets/*.yaml'), ('New designs (examples/)', 'examples/*.yaml'),
          ('Variant matrix (examples/matrix/)', 'examples/matrix/*.yaml'),
          ('Decoupling grid (examples/decoupling/)', 'examples/decoupling/*.yaml')]


def caption(name):
    """Return a build's title for the sheet and the path of its layout figure, which is None when the record is
    missing or has no validation."""
    p = os.path.join(ROOT, 'runs', name, name + '.json')
    if not os.path.exists(p):
        return 'no record', None
    r = json.load(open(p))
    a = r['assessment']
    v = a['verdict']
    img = os.path.join(ROOT, 'runs', name, name + '_layout.png')
    if 'validation' not in r:
        return '%s\n%s' % (name, textwrap.shorten(v, 240)), None
    V = r['validation']
    cd, pc = V.get('components_drawn') or {}, V.get('part_clearance') or {}
    m = r['input'].get('meander', {})
    head = ('%s head: %d units, %d arms%s, gap 1 %s' % (
        m.get('style', 'pack'), cd.get('head_units', 0), cd.get('head_arms', 0),
        (' (lane pairs %s)' % cd['arm_lane_pairs']) if cd.get('arm_lane_pairs') else '',
        'sealed' if cd.get('gap1_sealed') else 'at %.2f f* (declared)' % (cd.get('gap1_f_over_fstar') or 0))
        if m.get('enabled') else 'no head (%s start)' % r['input'].get('start', {}).get('style', 'sealed'))
    comp = '%s; inlet %s (%d filter posts); outlet %s; %d leg(s)' % (
        head, r['input'].get('inlet', {}).get('style', 'plain'), cd.get('filter_posts_counted', 0),
        cd.get('outlet_style', '?'), cd.get('legs', 0))
    chk = {c['check']: c['ok'] for c in a.get('checks', [])}
    ni = chk.get('no intersections: parts keep the floor apart')
    cr = chk.get('drawn components as requested')
    line3 = 'no intersections: %s (closest unrelated parts %s um); components as requested: %s' % (
        {True: 'yes', False: 'NO', None: 'n/a'}[ni], (pc.get('closest') or [{}])[0].get('gap_um', '-'),
        {True: 'yes', False: 'NO', None: 'n/a'}[cr])
    return '%s  [exit %s]\n%s\n%s\n%s' % (name, a.get('exit_code'), textwrap.shorten(v, 200), comp, line3), img


def main(out):
    """Write the sheet, four builds to a page, group by group."""
    with PdfPages(out) as pdf:
        for title, pat in GROUPS:
            names = sorted(os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(ROOT, pat)))
            for i in range(0, len(names), 4):
                fig, axs = plt.subplots(4, 1, figsize=(11, 15))
                fig.suptitle(title, fontsize=12)
                for ax, n in zip(axs, names[i:i + 4] + [None] * 4):
                    ax.axis('off')
                    if n is None:
                        continue
                    cap, img = caption(n)
                    if img and os.path.exists(img):
                        ax.imshow(plt.imread(img))
                    ax.set_title(cap, fontsize=6.5, loc='left')
                fig.tight_layout(rect=(0, 0, 1, 0.98))
                pdf.savefig(fig, dpi=110)
                plt.close(fig)
    print(out)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'docs', 'CONTACT_SHEET.pdf'))
