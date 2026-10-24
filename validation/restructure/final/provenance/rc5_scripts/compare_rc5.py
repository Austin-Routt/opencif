"""Compare every rc5 record with its rc4 record (two runs/ directories): exit code, verdict, design, network, realised
cutoff, Q1, the outlet (lead width and length, the stretch-end term), the aspect-ratio check, the model band and the
outlet-balance scenario, the robustness scenarios, the mask's sha256.  Prints a block per record and a summary.
Usage: python3 compare_rc5.py RC4_RUNS RC5_RUNS"""
import glob, json, os, sys

A = sys.argv[1] if len(sys.argv) > 2 else 'opencif_rc4/runs'
B = sys.argv[2] if len(sys.argv) > 2 else 'opencif_rc5/runs'


def load(root, n):
    p = os.path.join(root, n, n + '.json')
    return json.load(open(p)) if os.path.exists(p) else None


def chk(a, prefix):
    for c in a['assessment']['checks']:
        if c['check'].startswith(prefix):
            return c
    return None


def scen(V, key):
    return (V.get('model_scenarios') or {}).get(key) or {}


names = sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(os.path.join(A, '*', '*.json'))
               if os.path.basename(os.path.dirname(p)) + '.json' == os.path.basename(p))
summary = dict(records=0, same_mask=[], changed_mask=[], exit_changes={}, verdict_changes=[], design_changes=[],
               network_changes=[], hashes=set(), missing=[], ar_moves=[])
for n in names:
    a, b = load(A, n), load(B, n)
    if b is None:
        summary['missing'].append(n)
        continue
    summary['records'] += 1
    summary['hashes'].add((b.get('code') or {}).get('sha256') or (b.get('element_code') or {}).get('sha256'))
    if 'element' in a:                                  # arrays
        print('%-48s ARRAY exit %s -> %s | %s' % (n, a['exit_code'], b['exit_code'], b['verdict'][:100]))
        if a['exit_code'] != b['exit_code']:
            summary['exit_changes'][n] = (a['exit_code'], b['exit_code'])
        continue
    aa, ab = a['assessment'], b['assessment']
    lines = []
    if aa['exit_code'] != ab['exit_code']:
        summary['exit_changes'][n] = (aa['exit_code'], ab['exit_code'])
    if aa['verdict'] != ab['verdict']:
        summary['verdict_changes'].append(n)
        lines.append('verdict rc4: ' + aa['verdict'][-170:])
        lines.append('verdict rc5: ' + ab['verdict'][-170:])
    if 'refused' in a or 'refused' in b:
        print('%-48s REFUSED exit %s -> %s %s' % (n, aa['exit_code'], ab['exit_code'],
                                                 'same reasons' if (a.get('refused') or {}).get('reasons') == (b.get('refused') or {}).get('reasons') else 'REASONS DIFFER'))
        for l in lines:
            print('      ', l)
        continue
    da, db = a['design'], b['design']
    dd = [k for k in ('N', 'NM', 'layout', 'n_legs', 'A', 'F') if da.get(k) != db.get(k)]
    if dd:
        summary['design_changes'].append((n, [(k, da.get(k), db.get(k)) for k in dd]))
        lines.append('design: ' + str([(k, da.get(k), db.get(k)) for k in dd]))
    Va, Vb = a['validation'], b['validation']
    ra, rb = Va['realised_cutoff'], Vb['realised_cutoff']
    netd = []
    for k in ('F_as_drawn', 'F_intended'):
        x, y = Va['network'].get(k), Vb['network'].get(k)
        if x is None or y is None or abs(x - y) > 1e-9:
            netd.append((k, x, y))
    for k in ('D_realised_min_um', 'D_realised_max_um', 'share_of_F_within_1%'):
        if abs(ra[k] - rb[k]) > 1e-9:
            netd.append((k, ra[k], rb[k]))
    for k in ('undeclared_rows_outside', 'declared_rows_outside', 'rows_outside_1%'):
        if ra.get(k) != rb.get(k):
            netd.append((k, ra.get(k), rb.get(k)))
    q1a, q1b = Va['drawing_fidelity']['cutoff_dev_max'], Vb['drawing_fidelity']['cutoff_dev_max']
    if abs(q1a - q1b) > 1e-9:
        netd.append(('Q1', q1a, q1b))
    if netd:
        summary['network_changes'].append(n)
        lines.append('network: ' + '; '.join('%s %s -> %s' % (k, ('%.6g' % x) if isinstance(x, float) else x,
                                                              ('%.6g' % y) if isinstance(y, float) else y) for k, x, y in netd))
    oa, ob = Va.get('outlet_measured') or {}, Vb.get('outlet_measured') or {}
    lines.append('outlet rel A rc4 %.2e rc5 %.2e; stretch_end %s; lead min %s' % (
        oa.get('RsA_out_rel', float('nan')), ob.get('RsA_out_rel', float('nan')),
        {k: round(v['c'], 4) for k, v in (ob.get('stretch_end') or {}).items()},
        {k: round(v, 2) for k, v in (ob.get('lead_width_min_um') or {}).items()}))
    ca, cb = chk(a, 'aspect ratio of the narrowest'), chk(b, 'aspect ratio of the narrowest')
    if ca and cb and (ca['value'] != cb['value'] or ca['ok'] != cb['ok']):
        lines.append('AR: rc4 %s %s -> rc5 %s %s' % (ca['ok'], ca['value'], cb['ok'], cb['value']))
        if ca['value'].get('at') != cb['value'].get('at') or abs(ca['value'].get('depth_over_width', 0) - cb['value'].get('depth_over_width', 0)) > 0.005:
            summary['ar_moves'].append((n, ca['value'], cb['value']))
    ja, jb = scen(Va, 'junction_band') or scen(Va, 'band'), scen(Vb, 'junction_band') or scen(Vb, 'band')
    ba, bb = chk(a, '1-D model band'), chk(b, '1-D model band')
    if ba and bb and ba['value'] != bb['value']:
        lines.append('band rc4: %s' % str(ba['value'])[:230])
        lines.append('band rc5: %s' % str(bb['value'])[:230])
    pa, pb = scen(Va, 'outlet_side_resistance_+1%'), scen(Vb, 'outlet_side_resistance_+1%')
    ma, mb = scen(Va, 'outlet_side_resistance_-1%'), scen(Vb, 'outlet_side_resistance_-1%')
    if pa and pb:
        lines.append('outlet +-1%%: rows %s/%s -> %s/%s; max %.3f/%.3f -> %.3f/%.3f D; lost %s -> %s' % (
            pa.get('rows_newly_outside_tol'), ma.get('rows_newly_outside_tol'), pb.get('rows_newly_outside_tol'),
            mb.get('rows_newly_outside_tol'), pa.get('max_cutoff_change') or 0, ma.get('max_cutoff_change') or 0,
            pb.get('max_cutoff_change') or 0, mb.get('max_cutoff_change') or 0, pa.get('rows_losing_their_cutoff'),
            pb.get('rows_losing_their_cutoff')))
    if Va.get('robustness_design_network') != Vb.get('robustness_design_network'):
        lines.append('robustness_design_network differs')
    sa, sb = (a.get('files') or {}).get('dxf_sha256'), (b.get('files') or {}).get('dxf_sha256')
    (summary['same_mask'] if sa == sb and sa else summary['changed_mask']).append(n)
    print('%-48s exit %s -> %s  mask %s' % (n, aa['exit_code'], ab['exit_code'], 'SAME' if sa == sb else 'CHANGED'))
    for l in lines:
        print('      ', l)
print()
print('RECORDS', summary['records'], 'MISSING', summary['missing'])
print('HASHES', sorted(x for x in summary['hashes'] if x))
print('EXIT CHANGES', summary['exit_changes'])
print('MASKS SAME %d: %s' % (len(summary['same_mask']), ', '.join(summary['same_mask'])))
print('MASKS CHANGED %d: %s' % (len(summary['changed_mask']), ', '.join(summary['changed_mask'])))
print('VERDICT CHANGES %d: %s' % (len(summary['verdict_changes']), ', '.join(summary['verdict_changes'])))
print('DESIGN CHANGES', summary['design_changes'])
print('NETWORK CHANGES %d: %s' % (len(summary['network_changes']), ', '.join(summary['network_changes'])))
print('AR MOVES', summary['ar_moves'])
