"""rc7 (register INV-1; code review R2): the figure-aligned Dinh preset's head is the depicted preset's head exactly as
the planner makes it (its depicted rows from data/depicted.json['dinh'], then rows at the printed f* until the head
seals): this script plans presets/dinh2024_v1_depicted.yaml and writes the head's gap fractions and row count into
presets/dinh2024_v1_figure_aligned.yaml (reproduction.f_head, reproduction.NM_force), so that the two presets cannot
drift apart.  Planner only (no drawing).
usage: python3 -I provenance/rc7_scripts/make_dinh_figure_aligned_head.py [--check]   (--check: exit 1 if it would change)"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from opencif import spec as SP                              # noqa: E402
from opencif.physics import schedule as S  # noqa: E402

sp, _ = SP.to_engine(SP.load(os.path.join(ROOT, 'presets', 'dinh2024_v1_depicted.yaml')))
S.clear_caches()
sc = S.plan(sp)
NM = int(sc.NM)
f = [float(v) for v in sc.f[:NM]]
P = os.path.join(ROOT, 'presets', 'dinh2024_v1_figure_aligned.yaml')
with open(P) as fh:
    txt = fh.read()
new = re.sub(r'(?m)^  NM_force: \d+$', '  NM_force: %d' % NM, txt)
line = ('  f_head: [%s]   # the depicted preset\'s planned head (provenance/rc7_scripts/make_dinh_figure_aligned_head.py)'
        % ', '.join(repr(v) for v in f))
new, n_ = re.subn(r'(?m)^  f_head: \[[^\]\n]*\][^\n]*$', lambda m_: line, new)
assert n_ == 1 and new.count('  NM_force: %d' % NM) == 1 and line in new
print('depicted plan: head rows %d (%d depicted, %d at the printed f*), head %.2f mm; f/f* rows 1-18: %s' % (
    NM, 18, NM - 18, sc.x[NM] / 1000.0, ' '.join('%.2f' % (v / sc.info['fstar']) for v in f[:18])))
if new == txt:
    print('unchanged')
elif '--check' in sys.argv:
    print('WOULD CHANGE', P)
    sys.exit(1)
else:
    with open(P, 'w') as fh:
        fh.write(new)
    print('written', P)
