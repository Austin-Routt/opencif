"""collect.py: collect every run record into one table, written as Markdown and as CSV.

Usage: python3 tools/collect.py OUT_PREFIX [names...]

The table has one row for each run under runs/ (or for each name given). The call writes OUT_PREFIX.md and
OUT_PREFIX.csv and prints the Markdown."""
import csv, glob, json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = sys.argv[1]
names = sys.argv[2:] or sorted(os.path.basename(p) for p in glob.glob(os.path.join(ROOT, 'runs', '*')))
rows = []


def _bend(ms):
    """Summarize a record's junction band as the rows newly outside the tolerance at each end and the largest change."""
    a, b = ms.get('junction_band_lo'), ms.get('junction_band_hi')
    if not a or not b:
        return '-'
    return '%d/%d rows, %.1f%% D' % (a['rows_newly_outside_tol'], b['rows_newly_outside_tol'],
                                     100 * max(a['max_cutoff_change'] or 0, b['max_cutoff_change'] or 0))


for n in names:
    p = os.path.join(ROOT, 'runs', n, n + '.json')
    if not os.path.exists(p):
        continue
    r = json.load(open(p))
    a = r.get('assessment', {})
    if 'refused' in r:
        rows.append(dict(name=n, verdict=a.get('verdict', 'REFUSED'), N=r['design']['N'], D=round(r['design']['D_star_um'], 4)))
        continue
    V, d = r['validation'], r['design']
    rc, net, fid = V['realised_cutoff'], V['network'], V['drawing_fidelity']
    drc = V.get('drc', {})
    rows.append(dict(
        name=n, verdict=a.get('verdict', '')[:90], layout='%s x%d' % (d['layout'], d['n_legs']), N=d['N'], NM=d['NM'],
        footprint_mm='%.1f x %.1f' % tuple(V['drawn_parameters']['footprint_mm']), D=round(rc['D_target_um'], 4),
        D_min=round(rc['D_realised_min_um'], 4), D_max=round(rc['D_realised_max_um'], 4),
        undeclared_outside=rc['undeclared_rows_outside'],
        declared='; '.join('%s rows %s' % (k.split(' (')[0], v[:2]) for k, v in rc['declared_rows'].items()) or '-',
        share_within_1pct=round(100 * rc['share_of_F_within_1%'], 2), F_target=net['F_target'], F_drawn=round(net['F_as_drawn'], 5),
        Q1_cutoff_dev=float('%.2g' % fid['cutoff_dev_max']), narrowest_um=round(V['side_continuity']['narrowest_um'], 2),
        drc='ok' if drc and drc.get('n_lost_open') == 0 and drc.get('n_lost_solid') == 0 and drc.get('n_components') == 1 else str(drc.get('error', drc))[:30],
        envelope='ok' if (V['gates']['plateau_ok'] and V['gates']['M13_ok'] and V['gates']['Re_ok']
                          and V['gates'].get('stokes_ok', False)) else
                 'u %.3f M13 %.3f Re %.1f Re_gap %s Re_delta %s' % (
                     V['gates']['plateau_u_min'], V['gates']['M13_gap_dominance_max'], V['gates']['Re_max'],
                     ('%.1e' % V['gates']['Re_gap_max']) if 'Re_gap_max' in V['gates'] else 'n/a (old record)',
                     ('%.1e' % V['gates']['Re_delta_max']) if 'Re_delta_max' in V['gates'] else 'n/a'),
        bend_band=_bend(V.get('model_scenarios', {})),
        code=r.get('code', {}).get('sha256', 'none'),
        minutes=round(r.get('seconds_total', r.get('seconds', 0)) / 60, 1)))
keys = list(dict.fromkeys(k for r in rows for k in r))
with open(out + '.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, keys); w.writeheader(); [w.writerow(r) for r in rows]
with open(out + '.md', 'w') as fh:
    fh.write('| ' + ' | '.join(keys) + ' |\n|' + '---|' * len(keys) + '\n')
    for r in rows:
        fh.write('| ' + ' | '.join(str(r.get(k, '')) for k in keys) + ' |\n')
print(open(out + '.md').read())
