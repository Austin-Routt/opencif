"""Apply the report's §12.1 documentation edits (the Dinh row pitch read at the bars' anti-aliased length, 45.1 µm,
and the length paragraph's numbers at that pitch) to docs/DISSERTATION_NOTES.md and docs/REPRODUCTIONS.md.
Whitespace-tolerant (the docs wrap their sentences); every edit must match exactly once; files are rewritten as new
files (safe in a hard-linked working copy).  Not hashed: the code stamp is unchanged.

usage: python3 apply_doc_edits.py REPO_ROOT"""
import os
import re
import sys

ROOT = sys.argv[1]

PITCH_NEW = ('45.1 ± 0.1 µm pitch measured on Fig. 2B\'s second, third and fourth panels (45.16, 45.09, 45.08 µm; each '
             'panel\'s 50 µm bar read with its anti-aliased ends, 46.91 px; rc5 read 46 px and 46.0 ± 1.0 µm)')
COMMON = [
    ('at least 1897 rows, 87.3 mm,', 'at least 1897 rows, 85.6 mm,'),
    ('76.1 mm at F 0.5105 and 74.9 mm at the measured F 0.505', '74.6 mm at F 0.5105 and 73.4 mm at the measured F 0.505'),
    ('107.9 mm of rows for a width linear in the row index', '105.8 mm of rows for a width linear in the row index'),
    ('127.6–130.5 mm of head and rows (101.9 mm of linear rows)', '125.6–128.5 mm of head and rows (99.9 mm of linear rows)'),
    ('at least 98.8 % of f_max, as a constant 52.3 µm centre would', 'at least 96.9 % of f_max, as a constant 52.8 µm centre would'),
]
EDITS = {
    'docs/DISSERTATION_NOTES.md': [
        ('every gap row at least the 46.0 ± 1.0 µm pitch measured on Fig. 2B.', 'every gap row at least the ' + PITCH_NEW + '.'),
        ('At the measured F that bound leaves 0.9 mm for the inlet, filter and outlets',
         'At the measured F that bound leaves 2.4 mm for the inlet, filter and outlets'),
        ('credited, the bound is 74.1–74.5 mm, leaving up to 1.7 mm', 'credited, the bound is 72.6–73.0 mm, leaving up to 3.2 mm'),
        ('(a linear narrowing then 52 µm: within the first 3.8 % of the rows; 10.9 % at a 45 µm pitch, none at 47 µm)',
         '(a linear narrowing then 52 µm: within the first 10.2 % of the rows, 13.6 % crediting the depicted head; 8.8–11.7 % '
         'over the pitch\'s ± 0.2 µm)'),
        ('a linear narrowing confined to the first third of the rows would fit',
         'a linear narrowing confined to the first 39 % of the rows would fit'),
        ('the every-row bound alone 134.7–145.1 mm', 'the every-row bound alone 132.1–142.3 mm'),
        ('(a linear taper still needs 93.4 mm at the shortest)', '(a linear taper still needs 91.6 mm at the shortest)'),
    ] + COMMON,
    'docs/REPRODUCTIONS.md': [
        ('46.0 ± 1.0 µm pitch measured on Fig. 2B', PITCH_NEW),
        ('74.1–74.5 mm', '72.6–73.0 mm'),
        ('only 0.9 mm (up to 1.7 mm)', 'only 2.4 mm (up to 3.2 mm)'),
        ('within the first 3.8 % of the rows at the measured pitch; 10.9 % at 45 µm, none at 47 µm',
         'within the first 10.2 % of the rows at the measured pitch; 8.8–11.7 % over its ± 0.2 µm'),
        ('first third of the rows would fit', 'first 39 % of the rows would fit'),
    ] + COMMON,
}


def pattern(s):
    return re.compile(r'\s+'.join(re.escape(w) for w in s.split()))


for rel, edits in EDITS.items():
    p = os.path.join(ROOT, rel)
    t = open(p, encoding='utf-8').read()
    for old, new in edits:
        rx = pattern(old)
        n = len(rx.findall(t))
        if n != 1:
            raise SystemExit('%s: %d matches for %r' % (rel, n, old[:60]))
        t = rx.sub(lambda m: new, t)
    os.remove(p)                       # a new inode: never write through a hard link
    open(p, 'w', encoding='utf-8').write(t)
    print('%s: %d edits' % (rel, len(edits)))
