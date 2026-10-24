"""rc6: the Dinh V1 presets' paper.note and paper.length_bound rewritten from the rc6 computations (every number is
read from the outputs, none typed): dinh_length_bound_{depicted,specified}.json, dinh_taper_alternatives.json and
dinh_schedule_at_F.json in provenance/rc6_scripts/ (the rc6 ports at the presets' own pitch G + g_p = 45.1 um).  The
rc5 wording is kept where its conclusion still holds; where a number changes a conclusion's sign (the every-row bound
at F 0.5105, 0.3 mm over the printed footprint at rc5's 46 um pitch, now under it) the wording follows the number.
tests/test_figure_anchors.py checks the block against the same outputs.  Usage: python3 apply_dinh_notes.py OPENCIF_ROOT"""
import json
import os
import sys
from decimal import Decimal, ROUND_HALF_UP

R = sys.argv[1]
H = os.path.join(R, 'provenance', 'rc6_scripts')
TAP = json.load(open(os.path.join(H, 'dinh_taper_alternatives.json')))
SCH = json.load(open(os.path.join(H, 'dinh_schedule_at_F.json')))
PRINTED = 75.8
tapF = {round(q['F'], 4): q for q in TAP['readings']}
sens = {(r['F'], r['pitch_um']): r for r in TAP['sensitivity']}
pitch, tol = TAP['pitch_um'], TAP['pitch_tolerance_um']
pl, pm, ph = round(pitch - tol, 6), round(pitch, 6), round(pitch + tol, 6)


def g(v):
    return '%g' % v


def r1(v):
    """one decimal, half up, of the value as the outputs print it (85.55 -> 85.6, not the binary 85.5)"""
    return str(Decimal(repr(round(v, 3))).quantize(Decimal('0.1'), rounding=ROUND_HALF_UP))


def pct(v):
    return '%.1f' % (100.0 * v)


for k in ('depicted', 'specified'):
    name = 'dinh2024_v1_' + k
    B = json.load(open(os.path.join(H, 'dinh_length_bound_%s.json' % k)))
    assert abs(B['pitch_um'] - pitch) < 1e-9 and abs(B['pitch_tolerance_um'] - tol) < 1e-9, 'bound and taper pitches differ'

    def bound(F, which):
        got = [b for b in B['bounds'] if abs(b['F'] - F) < 1e-4 and which in b['counted']]
        assert len(got) == 1
        return got[0]
    cD = {F: bound(F, 'constant D')['length_mm_at_pitch'] for F in (0.5593, 0.5105, 0.505)}
    oh = {F: bound(F, 'own head')['length_mm_at_pitch'] for F in (0.5593, 0.5105, 0.505)}
    N1 = bound(0.5593, 'constant D')['N_min']
    plans = {round(q['F_target'], 4): q for q in SCH[name]['plans']}
    hp = {F: plans[F]['head_plus_rows_mm'] for F in (0.5592, 0.5105, 0.505, 0.472)}
    lin = {F: tapF[F]['length_mm_linear_75_to_52'] for F in (0.5593, 0.5105, 0.505, 0.472)}
    meanf = tapF[0.505]['constant_f_needed_in_printed_length_over_f_max']
    cw = tapF[0.505]['the_same_as_a_constant_centre_width_um']
    s, sl, sh_ = (sens[(0.505, p)]['largest_share_tapering'] for p in (pm, pl, ph))
    sc = sens[(0.505, pm)]['largest_share_tapering_with_%s_head_credit' % k]
    smax = max(v for p in (pl, pm, ph) for v in (sens[(0.505, p)]['largest_share_tapering'],
                                                 sens[(0.505, p)]['largest_share_tapering_with_%s_head_credit' % k])
               if v is not None)
    s472 = sens[(0.472, pm)]['largest_share_tapering']
    assert s == tapF[0.505]['largest_share_of_rows_tapering_if_rows_had_the_whole_length'] and \
        s472 == tapF[0.472]['largest_share_of_rows_tapering_if_rows_had_the_whole_length']
    assert cD[0.5593] > PRINTED and cD[0.5105] < PRINTED and cD[0.505] < PRINTED and hp[0.505] > PRINTED and \
        lin[0.505] > PRINTED and s is not None and s472 is not None, 'the wording below assumes these signs'
    left = {F: PRINTED - cD[F] for F in cD}
    lefth = {F: PRINTED - oh[F] for F in oh}
    note = (
        "length (premises: the frozen cutoff law's cubic kernel with the lamina the same at every gap, at the plateau "
        "prefactor kappa_g, which is App. B's own reading of Dinh's rule ('Under the Xia-scaled reading, which is the one "
        "Dinh describes, the lamina is the same at every gap', App. B :2970-2972; each gap adjusted 'to maintain the "
        "desired size cut-off throughout the device', p. 916), fixed by the printed f* at 75 um, so kappa drops out of the "
        "row count; every centre width >= 52 um; two gaps per row; every gap row at least the %s +- %s um post pitch "
        "measured on Fig. 2B at its 50 um bars read with their anti-aliased ends, rc5 46.0 +- 1.0).  (1) A lower bound "
        "for any taper, every row at f_max: the gap rows alone need >= %d rows, %s mm, at the F this preset closes "
        "(0.5593, the 1-D outlet identity at the printed 44 / 52 um), %s mm at F 0.5105 and %s mm at the measured F "
        "0.505 (%s / %s / %s mm with this preset's head, whose first rows are drawn above f_max, a departure from "
        "constant D), against the printed 75.8 mm of 'An assembled CIF device' (Fig. 3A), which also holds the inlet with "
        "its filter and the two outlets (p. 916): the bound leaves them %s mm (%s with the head's credit) at the "
        "measured F and %s mm (%s) at 0.5105, so by itself it shows no conflict at either.  (2) A centre channel "
        "narrowed 'progressively' (p. 916) needs more rows: this preset's own schedule needs %s mm of head and rows at "
        "F 0.505 (%s mm of linear rows), a width linear in the row index %s mm of rows; the rows fit 75.8 mm only if "
        "they extract on average at least %s %% of f_max, as a constant %s um centre would (a linear narrowing then "
        "52 um: within the first %s %% of the rows at the measured pitch, %s %% at %s um, %s %% at %s um, %s %% with this "
        "preset's head credit).  So under these premises the conflict holds at this preset's F, and at the measured mean "
        "F for a centre narrowed over more than the first %s %% of the rows (at most %s %% within the pitch's tolerance "
        "and with the head's credit), this preset's schedule and a width linear in the row index included (at one SD "
        "below the measured mean, F 0.472, a linear narrowing within the first %s %% of the rows would fit); it is "
        "between the premised model and the printed footprint, not among the printed numbers, and its cause (the "
        "kernel's lamina-flux relation in 140 um deep channels, the constant lamina, the pitch, or a taper or design "
        "rule the paper does not print) is not known.  The cutoff law below its plateau (rows wider than 62 um) does not "
        "remove it: holding D itself constant with the Hill curve there either lengthens the device or does not return "
        "the printed f* (provenance/rc6_scripts/dinh_length_bound.py, dinh_schedule_at_F.py, "
        "dinh_taper_alternatives.py; rc5's at 46.0 um in provenance/rc5_scripts/)"
        % (g(pitch), g(tol), N1, r1(cD[0.5593]), r1(cD[0.5105]), r1(cD[0.505]), r1(oh[0.5593]), r1(oh[0.5105]),
           r1(oh[0.505]), r1(left[0.505]), r1(lefth[0.505]), r1(left[0.5105]), r1(lefth[0.5105]), r1(hp[0.505]),
           r1(plans[0.505]['linear_rows_mm']), r1(lin[0.505]), r1(100.0 * meanf), r1(cw), pct(s), pct(sl), g(pl),
           pct(sh_), g(ph), pct(sc), pct(s), pct(smax), pct(s472)))
    status = (
        "conditional on the premises: the every-row bound conflicts at the outlet identity's F (%s mm against 75.8) "
        "and leaves %s mm (%s with this head's credit) at 0.5105 and %s mm (%s) at the measured F for the inlet, filter "
        "and outlets; at the measured F and pitch a linear narrowing fits only within the first %s %% of the rows "
        "(%s %% with this head's credit), so a progressive taper, this preset's included, conflicts; at one SD below "
        "the measured F a linear narrowing within the first %s %% of the rows fits"
        % (r1(cD[0.5593]), r1(left[0.5105]), r1(lefth[0.5105]), r1(left[0.505]), r1(lefth[0.505]), pct(s), pct(sc),
           pct(s472)))

    def wrap(text, first, cont, width=120):
        """wrap at spaces, never inside 'Fig. 2B', 'p. 916', 'App. B' or before a '%'"""
        for a_ in ('Fig. ', 'p. ', 'App. '):
            text = text.replace(a_, a_[:-1] + '\x00')
        text = text.replace(' %', '\x00%')
        out, line = [], first
        for w in text.split(' '):
            if len(line) + 1 + len(w) > width and line.strip():
                out.append(line)
                line = cont + w
            else:
                line = (line + ' ' + w) if line not in (first, cont) else line + w
        out.append(line)
        return '\n'.join(out).replace('\x00', ' ')

    def m(d, fmt='F_%s'):
        return '{' + ', '.join('%s: %.2f' % (fmt % g(F), v) for F, v in d.items()) + '}'
    block = (wrap(note, '  note: "', '    ') + '"\n' +
             wrap(status, '  length_bound: {status: "', '      ') + '", pitch_um: %s,\n' % g(pitch) +
             '    pitch_tolerance_um: %s, f_max: %.4e, printed_footprint_mm: %s,\n' % (g(tol), B['f_max'], g(PRINTED)) +
             '    min_length_mm_constant_D_head: %s,\n' % m(cD) +
             '    min_length_mm_this_head: %s,\n' % m(oh) +
             '    this_schedule_head_plus_rows_mm: %s,\n' % m(hp) +
             '    linear_taper_rows_mm: %s,\n' % m(lin) +
             '    mean_f_needed_over_f_max_at_F_0.505: %.3f,\n' % meanf +
             '    largest_share_tapering_at_F_0.505: {pitch_%s_um: %s, pitch_%s_um: %s, pitch_%s_um: %s,\n'
             % (g(pl), sl, g(pm), s, g(ph), sh_) +
             '      pitch_%s_um_with_this_head_credit: %s},\n' % (g(pm), sc) +
             '    script: provenance/rc6_scripts/dinh_length_bound.py,\n'
             '    scripts_after_rc5_reviews: [provenance/rc6_scripts/dinh_schedule_at_F.py, '
             'provenance/rc6_scripts/dinh_taper_alternatives.py]}\n')
    block = block.replace('f_max: %.4e' % B['f_max'], 'f_max: %s' % ('%.4e' % B['f_max']).replace('e-04', 'e-4'))
    p = os.path.join(R, 'presets', name + '.yaml')
    t = open(p).read()
    a = t.index('  note: "length (premises')
    b = t.index('# flow: Dinh 2024 runs V1')
    t = t[:a] + block + t[b:]
    open(p, 'w').write(t)
    import yaml
    lb = yaml.safe_load(open(p))['paper']['length_bound']
    assert lb['min_length_mm_constant_D_head']['F_0.505'] == round(cD[0.505], 2)
    print(name, 'note and length_bound rewritten; F 0.505 bound %.2f mm, schedule %.2f mm, share %s' % (cD[0.505], hp[0.505], s))
