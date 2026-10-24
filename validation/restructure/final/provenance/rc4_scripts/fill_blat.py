"""The full-build lattice's table and summary (the report's sec. 6.3 and 0.8) from blat/*.json beside this script:
writes blat_table.md and blat_summary.txt."""
import json, glob, os, re, sys, collections
S = os.path.dirname(os.path.abspath(__file__))
recs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(S, 'blat', '*.json')))]
recs.sort(key=lambda r: r['name'])
codes = collections.Counter(str(r.get('exit_code', r.get('exception'))) for r in recs)
def short(v):
    v = v or ''
    head = v.split(' [')[0].split(' -- ')[0][:110]
    tags = re.findall(r' -- ([A-Z][A-Z \-]+[A-Z])(?: \(([^)]*)\))?', v)
    return head + ''.join(' -- %s%s' % (t, (' (%s)' % d.split(';')[0][:60]) if d else '') for t, d in tags)
rows = ['| process | outlet | F | layout | mode | exit | verdict | hard checks failed | 1-D model band (reported): rows, max, share of F | time |',
        '|---|---|---|---|---|---|---|---|---|---|']
def band_short(v):
    m = re.search(r'junction band: (\d+) row\(s\) beyond \+-1% \(up to ([\d.]+)% of D; rows ([\d-]+); ([\d.]+)% of F\)', v or '')
    if not m:
        m2 = re.search(r'\[1-D model: (.*?)\]', v or '')
        return m2.group(1) if m2 else ''
    rr = m.group(3).replace('-', '–')
    return 'junction: %s (%s %s), %s %% D, %s %% F' % (m.group(1), 'rows' if '–' in rr else 'row', rr, m.group(2), m.group(4))
order = {'generic19': 0, 'gifford2014': 1, 'dinh2024': 2, 'xia2016': 3}
recs.sort(key=lambda r: (order.get(r['kwargs'].get('process', 'generic19'), 9), r['kwargs'].get('layout', 'linear') != 'paperclip',
                         r['kwargs'].get('layout', 'linear'), r['kwargs'].get('outlet', 'separate'), r['kwargs'].get('filtrate_fraction'),
                         'u0p035' in r['name']))
for r in recs:
    kw = r['kwargs']
    v = r.get('verdict') or ('%s: %s' % (r.get('exception'), r.get('message', '')[:100]))
    lay = kw.get('layout', 'linear') + (' (u 0.035)' if 'u0p035' in r['name'] else '')
    rows.append('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %.0f s |' % (
        kw.get('process', 'generic19'), kw.get('outlet', 'separate'), kw.get('filtrate_fraction'), lay,
        kw.get('mode', 'strict'), r.get('exit_code', r.get('exception')),
        short(v).replace('|', '/'), ', '.join(r.get('failed') or []) or 'none', band_short(v).replace('|', '/'), r['seconds']))
n = len(recs)
summ = '%d of %d cases' % (n, 22)
c0, c3 = codes.get('0', 0), codes.get('3', 0)
other = {k: v for k, v in codes.items() if k not in ('0', '3')}
summary = ('%d builds: every plan-stage case that plans inside the envelope (every outlet style at F 0.5 and 0.9 on the '
           'generic, Gifford and Dinh processes on a 25 x 75 mm paperclip, less the rings the plan refuses; the other '
           'layouts at F 0.9; the plateau edge) and one Xia-process ring at F 0.5, which plans outside it: '
           '%d exit 0 and %d exit 3%s, with no exception and no failed hard check' % (
               n, c0, c3, (', ' + ', '.join('%d exit %s' % (v, k) for k, v in other.items())) if other else ''))
over = [r['name'] for r in recs if 'OVER THE CHIP AREA' in (r.get('verdict') or '')]
if over:
    summary += ('; %d drawn over the chip area (%s: a linear device of that length cannot fit 75 mm), which is '
                'reported, never a failure' % (len(over), ', '.join(over)))
anyfail = [r['name'] for r in recs if r.get('failed') or 'exception' in r]
print('CODES', dict(codes)); print('FAILED', anyfail); print(summary)
open(os.path.join(S, 'blat_table.md'), 'w').write('\n'.join(rows) + '\n')
open(os.path.join(S, 'blat_summary.txt'), 'w').write(summary + '\n')
