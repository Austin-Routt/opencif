"""results.py: write docs/RESULTS.md and docs/results.csv from the run records, with one row for each preset or
example, grouped as in GROUPS.

Usage: python3 tools/results.py   (reads runs/<name>/<name>.json for every spec under presets/ and examples/)"""
import csv
import glob
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GROUPS = [('Reproductions (presets/)', 'presets/*.yaml'), ('New designs (examples/)', 'examples/*.yaml'),
          ('Variant matrix (examples/matrix/)', 'examples/matrix/*.yaml'),
          ('Decoupling grid (examples/decoupling/)', 'examples/decoupling/*.yaml')]


# the specs built to show a check catching a real departure, whose FAIL or REFUSED is the expected outcome
EXPECTED = {'m14_paperclip_const_nomeander_floorstart19': 'FAIL by design, because a floor start built as a new design '
            'does not exempt its over-extracting entrance rows',
            'wbc_xia_process_refused': 'REFUSED by design, because the Xia process (M13 0.18) is outside the envelope in strict mode',
            'dinh_targets_sealed_head': 'exit 3 by design, for Dinh\'s targets at Dinh\'s width (u 0.029, below the plateau), built in '
            'evaluate mode to show the sealed head and the C7-method band'}


def rows_txt(v):
    """Return a count of rows with their first and last row, as '<count> (<first>-<last>)', or '0'."""
    return '%d (%d-%d)' % (v[2], v[0], v[1]) if v else '0'


def load(name):
    """Return the run record of a spec, or None when it has none."""
    p = os.path.join(ROOT, 'runs', name, name + '.json')
    return json.load(open(p)) if os.path.exists(p) else None


def env_txt(g):
    """Return the envelope's numbers in compact form, which are u_min, the gap dominance, the channel Re, and the larger
    Stokes number."""
    if not g:
        return '-'
    st = max(g.get('Re_gap_max', 0.0), g.get('Re_delta_max', 0.0))
    return 'u %.4f; M13 %.3f; Re %.1f; Re_st %.1e' % (g['plateau_u_min'], g['M13_gap_dominance_max'], g['Re_max'], st)


def row(name):
    """Return the table row of one spec from its record."""
    r = load(name)
    if r is None:
        return dict(name=name, verdict='(no record, because the spec was not built or the generator raised; see logs/%s.log)' % name, code='-')
    a = r.get('assessment', {})
    out = dict(name=name, verdict=a.get('verdict', ''), exit=a.get('exit_code'), code=r.get('code', {}).get('sha256', 'none'))
    if name in EXPECTED:
        out['verdict'] = '%s (expected: %s)' % (out['verdict'], EXPECTED[name])
    d = r['design']
    if 'refused' in r:
        g = r['refused']['gates']
        out.update(N=d['N'], NM=d.get('NM'), D='%.4f' % d['D_star_um'], F=('%.5f' % d['F']) if d.get('F') is not None else '-',
                   band=(('%.3f-%.3f' % tuple(g['cutoff_band_um']['band_um'])) if g['cutoff_band_um'].get('band_um') else
                         'none') if g else '-', envelope=env_txt(g))
        return out
    V = r['validation']
    rc, net, fid = V['realised_cutoff'], V['network'], V['drawing_fidelity']
    dec = rc.get('declared_by_option') or {}
    out.update(layout='%s x%d' % (d['layout'], d['n_legs']), N=d['N'], NM=d['NM'],
               footprint='%.1f x %.1f' % tuple(V['drawn_parameters']['footprint_mm']), D='%.4f' % rc['D_target_um'],
               undeclared_out=rc['undeclared_rows_outside'],
               within1='%.2f' % (100 * rc['share_of_F_within_1%']),
               F='%.5f / %.5f' % (net['F_target'], net['F_as_drawn']), Q1='%.1e' % fid['cutoff_dev_max'],
               declared='; '.join('%s: %d rows%s' % (k, v['rows_outside'], (', worst %.2f D' % v['worst_over']) if v.get('worst_over') else '')
                                  for k, v in dec.items() if v['rows_outside']) or '-',
               narrowest='%.2f' % V['side_continuity']['narrowest_um'], envelope=env_txt(V.get('gates')),
               aspect='%.1f' % (r['input']['process']['depth_um'] / V['side_continuity']['narrowest_um']))
    # the aspect-ratio check's own value (the narrowest measured CIF channel of any type), where the record has it
    for c_ in r['assessment'].get('checks', []):
        if c_['check'].startswith('aspect ratio of the narrowest') and isinstance(c_.get('value'), dict):
            out['aspect'] = '%.1f (%s)' % (c_['value']['depth_over_width'], c_['value']['at'])
    cd = V.get('components_drawn') or {}
    m = r['input'].get('meander', {})
    asked = m.get('style', 'pack') if m.get('enabled') else 'none'
    out['head'] = ('%s: %s%s' % (asked, cd.get('head_drawn', '-'),
                                 (', gap 1 sealed' if cd.get('gap1_sealed') else ', gap 1 declared') if cd.get('head_units') else '')
                   if cd else '-')
    chk = {c['check']: c['ok'] for c in a.get('checks', [])}
    out['parts'] = {True: 'ok', False: 'FAIL', None: '-'}[chk.get('no intersections: parts keep the floor apart')]
    out['as_asked'] = {True: 'ok', False: 'FAIL', None: '-'}[chk.get('drawn components as requested')]
    rb = V.get('robustness_design_network') or {}
    if all('rows_outside_1pct' in (rb.get(k) or {}) for k in ('bias +0.5 um', 'bias -0.5 um')):
        bp, bm = rb['bias +0.5 um'], rb['bias -0.5 um']
        lo_, hi_ = sorted((bp['rows_outside_1pct'], bm['rows_outside_1pct']))
        ext = [x for x in (bp.get('cutoff_row1_over_Dstar'), bm.get('cutoff_row1_over_Dstar'), bp['cutoff_max_over_Dstar'],
                           bm['cutoff_max_over_Dstar'], bp['cutoff_min_over_Dstar'], bm['cutoff_min_over_Dstar']) if x is not None]
        nl_ = max(bp.get('rows_without_cutoff', 0), bm.get('rows_without_cutoff', 0))
        out['bias'] = '%d-%d (%.1f%%); %.2f-%.2f D%s' % (lo_, hi_, 100.0 * hi_ / d['N'], min(ext), max(ext),
                                                       ('; %d lose their cutoff' % nl_) if nl_ else '')
    nz = rb.get('noise 0.010 um')
    if nz and 'rows_outside_1pct' in nz:
        out['noise'] = '%d' % nz['rows_outside_1pct']
    ms = V.get('model_scenarios') or {}
    b = [ms.get(k) for k in ('junction_band_lo', 'junction_band_hi') if ms.get(k)]
    if b:
        nl = max(x.get('rows_losing_their_cutoff', 0) for x in b)
        out['junction_band'] = '%d rows%s, %.1f%% D' % (max(x['rows_newly_outside_tol'] for x in b),
                                                        (' (%d lose cutoff)' % nl) if nl else '',
                                                        100 * max((x['max_cutoff_change'] or 0) for x in b))
    ob = [ms.get(k) for k in ('outlet_side_resistance_+1%', 'outlet_side_resistance_-1%') if ms.get(k)]
    if ob:
        out['outlet_balance'] = '%d rows, %.0f%% D, %.2f%% F' % (max(x['rows_newly_outside_tol'] for x in ob),
                                                               100 * max((x['max_cutoff_change'] or 0) for x in ob),
                                                               100 * max((x.get('share_of_F_in_those_rows') or 0) for x in ob))
    return out


def main():
    """Write docs/RESULTS.md and docs/results.csv, and print the Markdown."""
    lines = ['# Results', '',
             'Every record below was produced by `python3 -m opencif build <spec>`, and `code` is the hash of the '
             'package source and data that built it, which all rows should share. The columns are the exit code (0 '
             'PASS/REPRODUCED, 2 FAIL, 3 outside the validity envelope, 4 INCOMPLETE, 5 REFUSED); N gaps (NM in the '
             'head); the target D; the undeclared gaps outside D ± 1%; the share of the net filtrate extracted within '
             'D ± 1%; F, target and drawn; Q1 (the largest difference between the drawn and the intended cutoff, in '
             'units of D); the declared departures for each option; the head asked for and what was drawn (the arms '
             'and their lane pairs, the unit widths, and gap 1 sealed or declared); the no-intersection check (every '
             'part at least the floor from every part that it is not joined to) and the check that the drawn '
             'components are those requested; the junction model band (the rows that leave D ± 1% when the bend, head '
             'T-junction, and outlet stretch-end terms take either end of their range, and the largest change in units '
             'of D); the validity-envelope numbers (u_min = min D/w_c against the plateau 0.035, the gap dominance M13 '
             'against 0.05, the channel Re against 10, and the larger of the Stokes numbers Re_gap and Re_delta '
             "against 1e-2); the outlet balance (the rows that leave D ± 1% when both side outlets' resistance moves "
             "by ±1%, the largest change in units of D, and those rows' share of F, docs/LIMITATIONS.md item 15a); the "
             'depth over the narrowest drawn CIF channel (measured, since rc5, over gap openings, the center channel, '
             'side channels, U-turn sides, head units, and outlet leads, which is a gate in strict mode and is '
             'reported in the other modes, docs/LIMITATIONS.md item 7); and the fabrication sensitivity on the design '
             'network (validation.robustness_design_network). The bias and noise columns give the gaps outside D ± 1% '
             'under a uniform width bias of +0.5 and -0.5 µm on every lateral width (and their share of N), with the '
             'range of the cutoff over D across both biases (row 1 and the inner rows) and the gaps that lose their '
             'cutoff (reverse flow or f >= 1/2; a seam row), and the gaps outside D ± 1% for one draw (seed 1) of 10 '
             'nm random per-segment side-width noise. The ± 1% columns are statements about the nominal geometry, '
             'whereas these two say how far a real chip departs from it (docs/LIMITATIONS.md items 13-14).',
             '']
    allrows = []
    for title, pat in GROUPS:
        names = sorted(os.path.splitext(os.path.basename(p))[0] for p in glob.glob(os.path.join(ROOT, pat)))
        rs = [row(n) for n in names]
        allrows += rs
        lines += ['## ' + title, '', '| spec | exit | verdict | head drawn | no intersections | as asked | N (NM) | D µm | undeclared out | % F within 1 % | F target / drawn | Q1 | declared | junction band | outlet ±1 % | envelope | h / narrowest CIF channel | ±0.5 µm bias: gaps out (of N); cutoff/D | 10 nm noise: gaps out | code |',
                  '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
        for r in rs:
            lines.append('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
                r['name'], r.get('exit', '-'), r['verdict'].replace('|', '/'), r.get('head', '-'), r.get('parts', '-'),
                r.get('as_asked', '-'),
                ('%s (%s)' % (r['N'], r.get('NM', 0))) if 'N' in r else '-', r.get('D', '-'), r.get('undeclared_out', '-'),
                r.get('within1', '-'), r.get('F', '-'), r.get('Q1', '-'), r.get('declared', '-'), r.get('junction_band', '-'),
                r.get('outlet_balance', '-'), r.get('envelope', '-'), r.get('aspect', '-'), r.get('bias', '-'), r.get('noise', '-'), r['code']))
        lines.append('')
    open(os.path.join(ROOT, 'docs', 'RESULTS.md'), 'w').write('\n'.join(lines) + '\n')
    keys = list(dict.fromkeys(k for r in allrows for k in r))
    with open(os.path.join(ROOT, 'docs', 'results.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        for r in allrows:
            w.writerow(r)
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
