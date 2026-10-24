"""test_figure_anchors.py: the presets, and the drawn masks where they are built, against the papers' primary figures.

The suite answers the reviews of rc4. It compares the presets and the drawn masks with local measurements made on the
papers' own figures and calibrated on their scale bars (opencif/data/r6/figure_anchors.json). The rc6 reading measures
every bar with its anti-aliased ends (provenance/rc6_scripts/measure_rc6.py), the rc5 readings are kept in the file's
rc5_values block, and opencif/data/r5/figure_anchors.json is unchanged. This is not the preset/factory round trip of
test_presets_match_replicas.py, and it does not compare global schematic panels, which carry no scale (the anchors
file lists them as not_to_scale).

Against Dinh 2024 Fig. 2B (the outlet panel), both V1 presets draw the flare at the measured angle and straight run,
their post pitch G + g_p is the measured pitch, and the printed outlet widths, depth, and gap are the preset's. In the
drawn Dinh V1 masks (runs/, when built), the last side channel is the printed 44 um. The outer wall of each side
channel, measured on the fluid polygon, leaves the side channel's width at flare_start_um from the last gap's center
(the virtual corner, up to which the wall is within 0.25 um of straight), diverges at flare_deg, and follows the
designed profile to the stretch's end, and the record's design.outlet.L_out_um is that drawn stretch. In every built
record (runs/), design.outlet.L_out_um is the layout's drawn outlet stretch. The presets' pitch (+-0.2 um), l_g,
center width, and final side width are compared with those of the Gifford 2014 Fig. 4 inset. The pitch and l_g of
both Fig. 6 presets lie within the ranges measured on the Gifford 2014 Fig. 6 micrograph at its printed 125 um center
(provenance/rc6_scripts/gifford_fig6_measurement.json, measure_gifford_fig6.py), and the measurement's raster is the
recorded one. Every global panel that the reproductions rest on is marked not_to_scale or scale-uncertain. Finally,
the Dinh presets' paper.length_bound numbers are the provenance computations they quote, as the rc5 reviews asked,
computed by provenance/rc6_scripts/ at the presets' own pitch G + g_p, which must be the measured pitch.

Run it with `python3 tests/test_figure_anchors.py` (about a minute, because it loads every run pickle). The exit codes
follow the rc5 reviews, which found the suite printing ALL OK on the evaluation archive while 17 of its 33 checks could
not run there:
  0   ALL OK, and every check ran
  1   FAILED
  4   INCOMPLETE, because a check could not run, either because runs/ lacks the run pickles or records it reads (the
      evaluation archive ships the records without the pickles) or because a provenance output that a check reads is
      missing
The suite takes one option:
  --source-anchors-only   runs only the checks that need no build (the presets against the figures, and the Dinh
                          presets' length notes against their computations) and exits 0 if they all ran and passed.
                          The drawn-mask and record checks are then not run, by request, even where runs/ has them,
                          and a check that could not run for another reason still makes the run INCOMPLETE (exit 4)."""
import argparse
import json
import math
import os
import pickle
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT]
from opencif import spec as SP                        # noqa: E402
from opencif import pipeline as _PL_CODE                # noqa: E402  (every suite's log names the code it ran on)
print('opencif %s (code %s)' % (_PL_CODE.CODE['version'], _PL_CODE.CODE['sha256']))
import numpy as np                                     # noqa: E402

_ap = argparse.ArgumentParser(description='presets and drawn masks against scale-bar measurements on the papers\' figures')
_ap.add_argument('--source-anchors-only', action='store_true',
                 help='run the presets-against-figures checks only; the drawn-mask and record checks are not run')
ARGS = _ap.parse_args()
A = json.load(open(os.path.join(ROOT, 'opencif', 'data', 'r6', 'figure_anchors.json')))
fails = []
n_run = [0]
not_run = []                       # (what, number of checks, why, by request)


def check(name, cond, info=''):
    n_run[0] += 1
    print('%-96s %s %s' % (name, 'ok ' if cond else 'FAIL', info))
    if not cond:
        fails.append(name)


def skip(what, n, why, by_request=False):
    """Record and print n checks that did not run, with the reason why. by_request marks checks left out on request
    (--source-anchors-only), which do not make the run INCOMPLETE."""
    not_run.append((what, n, why, by_request))
    print('  NOT RUN (%d check%s): %s: %s' % (n, '' if n == 1 else 's', what, why))


# Each drawn Dinh mask has 8 checks, which are w_s(N), the record's stretch, and the run, angle, and profile of each
# side.
DINH_MASK_CHECKS = 8


dn = A['dinh2024_fig2b_panel4']
lo, hi = dn['flare_deg_per_side']['range']
on = dn['flare_onset_from_last_gap_centre_um']
on_mid = 0.5 * (on['left'] + on['right'])
for name in ('dinh2024_v1_depicted', 'dinh2024_v1_specified'):
    sp, op = SP.to_engine(SP.load(os.path.join(ROOT, 'presets', name + '.yaml')))
    check('%s: flare_deg %.2f within the measured %.2f-%.2f deg (+-0.3)' % (name, sp.flare_deg, lo, hi),
          lo - 0.3 <= sp.flare_deg <= hi + 0.3)
    check('%s: flare_start %.1f um within the measured %.1f +- %.1f um' % (name, sp.flare_start, on_mid, on['tolerance']),
          abs(sp.flare_start - on_mid) <= on['tolerance'])
    pitch = sp.G + sp.g_p
    check('%s: post pitch G + g_p = %.1f um within the measured %.2f +- %.1f um' % (
        name, pitch, dn['post_pitch_um']['value'], dn['post_pitch_um']['tolerance']),
        abs(pitch - dn['post_pitch_um']['value']) <= dn['post_pitch_um']['tolerance'])
    check('%s: the printed outlet widths (w_c %g, depth %g, G %g) are the preset\'s (%g, %g, %g)' % (
          name, dn['printed']['w_c_out_um'], dn['printed']['depth_um'], dn['printed']['G_um'], sp.wc_min, sp.h, sp.G),
          sp.wc_min == dn['printed']['w_c_out_um'] and sp.h == dn['printed']['depth_um'] and sp.G == dn['printed']['G_um'])

    # The drawn mask, when it is built, is measured along the outer wall of each side channel over the outlet stretch.
    pk = os.path.join(ROOT, 'runs', name, name + '.pkl')
    if ARGS.source_anchors_only:
        skip('%s drawn-mask checks' % name, DINH_MASK_CHECKS, 'not requested (--source-anchors-only)', True)
        continue
    if not os.path.exists(pk):
        skip('%s drawn-mask checks' % name, DINH_MASK_CHECKS, 'no run pickle %s' % os.path.relpath(pk, ROOT))
        continue
    from shapely.geometry import LineString
    d = pickle.load(open(pk, 'rb'))
    a, fluid = d['arrays'], d['fluid']
    N = int(a['N'])
    leg = a['legs'][-1]
    gx0 = float(a['gx'][-1])
    x_last = float(a['x'][N - 1]) - float(leg['p0'] or 0.0)
    y_in = float(a['wc'][-1]) / 2.0 + float(a['spec'].l_g)
    w_last = float(a['ws'][N - 2])
    spec_d = a['spec']
    check('%s drawn: the last side channel w_s(N) %.2f um is the printed %g um (+-0.1)' % (
        name, w_last, dn['printed']['w_s_out_um']), abs(w_last - dn['printed']['w_s_out_um']) <= 0.1)
    if getattr(spec_d, 'flare_deg', None) is None:
        check('%s drawn: the built mask carries the flare (a record built before rc5, or without the preset\'s '
              'outlet.flare_deg, does not)' % name, False)
        continue
    # The record's drawn stretch length (design.outlet.L_out_um) must be the layout's, which ends where the last leg
    # ends.
    rec_j = json.load(open(os.path.join(ROOT, 'runs', name, name + '.json')))
    L_rec = float(rec_j['design'].get('outlet', {}).get('L_out_um', float('nan')))
    L_lay = float(leg['p1']) - float(a['x'][N - 1])
    Lf0 = (float(a['outlet']['w_lead']) - w_last) / math.tan(math.radians(spec_d.flare_deg))
    check('%s drawn: the record\'s outlet stretch %.2f um is the layout\'s (%.2f), the straight run, the flare (%.1f) '
          'and the end arc' % (name, L_rec, L_lay, Lf0),
          abs(L_rec - L_lay) < 1e-6 and 0.0 <= L_rec - (spec_d.flare_start + Lf0) <= 25.0)
    for side, sg in (('A', 1.0), ('B', -1.0)):
        ss, ws = [], []
        for s in np.arange(2.0, L_lay - 2.0, 1.0):
            xg = gx0 + (x_last + s) * leg['dirn']
            cut = fluid.intersection(LineString([(xg, float(leg['Y']) - 3000.0), (xg, float(leg['Y']) + 3000.0)]))
            segs = [g for g in getattr(cut, 'geoms', [cut]) if g.length > 0]
            for g in segs:
                y0, y1 = sorted([(g.coords[0][1] - float(leg['Y'])) * sg * leg['dirn'],
                                 (g.coords[-1][1] - float(leg['Y'])) * sg * leg['dirn']])
                if abs(y0 - y_in) < 1.0:
                    ss.append(s)
                    ws.append(y1 - y_in)
        ss, ws = np.array(ss), np.array(ws)
        straight = ss[np.abs(ws - w_last) < 0.25]
        run = float(straight.max()) if straight.size else 0.0
        # The linear flare lies between its two tangent arcs, whose tangent length is about 3 um, so the fit skips 10 um
        # at each end.
        Lf = (float(a['outlet']['w_lead']) - w_last) / math.tan(math.radians(spec_d.flare_deg))
        m = (ss > spec_d.flare_start + 10.0) & (ss < spec_d.flare_start + Lf - 10.0)
        ang = math.degrees(math.atan(np.polyfit(ss[m], ws[m], 1)[0])) if m.sum() > 5 else float('nan')
        check('%s drawn, side %s: straight %.1f um from the last gap centre (preset %.1f; measured %.1f +- %.1f)' % (
            name, side, run, spec_d.flare_start, on_mid, on['tolerance']),
            abs(run - spec_d.flare_start) <= 4.0 and abs(run - on_mid) <= on['tolerance'] + 4.0)
        check('%s drawn, side %s: flare %.2f deg (preset %.2f; measured %.2f-%.2f)' % (
            name, side, ang, spec_d.flare_deg, lo, hi), abs(ang - spec_d.flare_deg) <= 0.2 and lo - 0.4 <= ang <= hi + 0.4)
        # The wall must follow the designed profile (straight run, arc, flare, and arc to the lead width) to the
        # stretch's end.
        import copy
        from opencif.geometry import layout as LY
        sp2 = copy.copy(spec_d)
        sp2.L_out = L_lay
        ps_, pw_ = np.array(LY.flare_profile(sp2, w_last, float(a['outlet']['w_lead']), L_lay)).T
        mm = ss > L_lay - 40.0
        dev = float(np.max(np.abs(np.interp(ss[mm], ps_, pw_) - ws[mm]))) if mm.any() else float('nan')
        check('%s drawn, side %s: the wall follows the designed profile over the stretch\'s last 40 um to %.1f um '
              '(largest difference %.3f um)' % (name, side, float(ss.max()) if ss.size else float('nan'), dev),
              dev < 0.25 and ss.size > 0 and ss.max() >= L_lay - 3.0)

# In every built record, design.outlet.L_out_um must be the drawn outlet stretch, the layout's (the last leg ends
# there), which a ring's loops, a turn, or the flare set. The rc5 review found the field reading the spec's base value.
import glob as _glob
n_rec, bad_rec, no_pk = 0, [], []
for jp in ([] if ARGS.source_anchors_only else sorted(_glob.glob(os.path.join(ROOT, 'runs', '*', '*.json')))):
    nm = os.path.basename(os.path.dirname(jp))
    pk_ = os.path.join(ROOT, 'runs', nm, nm + '.pkl')
    if os.path.basename(jp) != nm + '.json':
        continue
    rj = json.load(open(jp))
    if 'outlet' not in (rj.get('design') or {}):
        continue
    if not os.path.exists(pk_):
        # A record without its pickle cannot be checked, so it is counted rather than ignored.
        no_pk.append(nm)
        continue
    ar_ = pickle.load(open(pk_, 'rb'))['arrays']
    L_lay_ = float(ar_['legs'][-1]['p1']) - float(ar_['x'][int(ar_['N']) - 1])
    n_rec += 1
    if abs(float(rj['design']['outlet']['L_out_um']) - L_lay_) > 1e-6:
        bad_rec.append('%s %.2f vs %.2f' % (nm, rj['design']['outlet']['L_out_um'], L_lay_))
if n_rec:
    check('every built record (%d): design.outlet.L_out_um is the layout\'s drawn stretch' % n_rec, not bad_rec,
          '; '.join(bad_rec[:4]))
    if no_pk:                  # some records have pickles and some do not, so the check covered only part of runs/
        skip('the record-against-layout check for %d record(s) without run pickles' % len(no_pk), 1,
             'no run pickle: %s' % ', '.join(no_pk[:4]))
elif ARGS.source_anchors_only:
    skip('the record-against-layout check over every built record', 1, 'not requested (--source-anchors-only)', True)
else:                          # the evaluation archive ships the records without their run pickles
    skip('the record-against-layout check over every built record', 1, 'no run pickles in runs/')

gf = A['gifford2014_fig4_inset']
for name in ('gifford2014_fig4_design', 'gifford2014_fig4_as_fabricated_20um', 'gifford2014_fig4_as_fabricated_8um'):
    sp, op = SP.to_engine(SP.load(os.path.join(ROOT, 'presets', name + '.yaml')))
    check('%s: pitch G + g_p = %.2f um within the inset\'s %.2f +- 0.2 um' % (name, sp.G + sp.g_p, gf['row_pitch_um']),
          abs(sp.G + sp.g_p - gf['row_pitch_um']) <= 0.2)
    check('%s: l_g %.1f um is the inset\'s %.1f um (post columns centre to centre minus the printed centre)' % (
        name, sp.l_g, gf['l_g_um']), abs(sp.l_g - gf['l_g_um']) <= 0.05)
    check('%s: centre width %.0f um (printed) within the inset\'s %.1f +- %.1f um' % (
        name, sp.wc0, gf['centre_width_um'], gf['width_pick_uncertainty_um']),
        abs(sp.wc0 - gf['centre_width_um']) <= gf['width_pick_uncertainty_um'])
rec = os.path.join(ROOT, 'runs', 'gifford2014_fig4_design', 'gifford2014_fig4_design.json')
if ARGS.source_anchors_only:
    skip('gifford2014_fig4_design drawn final side width', 1, 'not requested (--source-anchors-only)', True)
elif not os.path.exists(rec):
    skip('gifford2014_fig4_design drawn final side width', 1, 'no record %s' % os.path.relpath(rec, ROOT))
else:
    w_out = json.load(open(rec))['validation']['drawn_parameters']['ws_out_eq']
    check('gifford2014_fig4_design drawn: final side width %.1f um within the inset\'s %.1f +- %.1f um' % (
        w_out, gf['final_side_width_um'], gf['width_pick_uncertainty_um']),
        abs(w_out - gf['final_side_width_um']) <= gf['width_pick_uncertainty_um'])
# Gifford Fig. 6 is the micrograph, measured at its printed center width, whereas rc5 used Fig. 4's l_g and an
# unscripted pitch.
G6P = os.path.join(ROOT, 'provenance', 'rc6_scripts', 'gifford_fig6_measurement.json')
if not os.path.exists(G6P):
    skip('Gifford Fig. 6 presets against the micrograph', 3, 'provenance output missing: %s' % os.path.relpath(G6P, ROOT))
else:
    G6 = json.load(open(G6P))
    check('gifford Fig. 6 measurement: made from the project raster gif7-000.png (recorded sha256 %s...)' % G6['raster_sha256'][:16],
          G6['raster_sha256'].startswith('182826e7aaccf14c'))
    ras6 = os.path.join(os.environ.get('OPENCIF_PAPER_RASTERS', ''), 'gif7-000.png')
    if os.environ.get('OPENCIF_PAPER_RASTERS') and os.path.exists(ras6):
        import hashlib
        h6 = hashlib.sha256(open(ras6, 'rb').read()).hexdigest()
        check('gifford Fig. 6: the raster in OPENCIF_PAPER_RASTERS is that one (sha256 %s...)' % h6[:16], h6 == G6['raster_sha256'])
    else:                                   # optional; the rasters are the publishers' and not in the repository
        print('  (optional, not counted: set OPENCIF_PAPER_RASTERS to a folder with gif7-000.png to hash the raster itself)')
    for name in ('gifford2014_fig6_design', 'gifford2014_fig6_as_fabricated_15um'):
        sp, op = SP.to_engine(SP.load(os.path.join(ROOT, 'presets', name + '.yaml')))
        lo_, hi_ = G6['row_pitch_range_um']
        llo, lhi = G6['l_g_range_um']
        check('%s: pitch G + g_p = %.2f um within Fig. 6\'s %.2f-%.2f um and l_g %.1f within %.1f-%.1f um' % (
            name, sp.G + sp.g_p, lo_, hi_, sp.l_g, llo, lhi),
            lo_ - 0.05 <= sp.G + sp.g_p <= hi_ + 0.05 and llo - 0.05 <= sp.l_g <= lhi + 0.05)
nts = A['not_to_scale']
check('global panels marked: Gifford Figs. 4 (left) and 6 not to scale; Dinh Fig. 2A main panel scale not established',
      all(k in nts for k in ('gifford2014_fig4_left', 'gifford2014_fig6', 'dinh2024_fig2a_main')))

# The Dinh presets' length notes are checked against the computations they quote, because a review of rc5 asked for a
# test that the wording's numbers are the computed ones. The computations are provenance/rc6_scripts/
# dinh_length_bound_*.json, dinh_schedule_at_F.json, and dinh_taper_alternatives.json, made at the presets' own pitch
# G + g_p (rc5's, at 46 um, stay in rc5_scripts/). The square-turn head's recomputation
# (provenance/rc7_scripts/dinh_lengths/, the rc6 scripts run on this tree) is read instead when it exists.
PV = os.path.join(ROOT, 'provenance', 'rc7_scripts', 'dinh_lengths')
if not os.path.isdir(PV):
    PV = os.path.join(ROOT, 'provenance', 'rc6_scripts')
for name in ('dinh2024_v1_depicted', 'dinh2024_v1_specified'):
    srcs = [os.path.join(PV, 'dinh_length_bound_%s.json' % name.split('_')[-1]), os.path.join(PV, 'dinh_schedule_at_F.json'),
            os.path.join(PV, 'dinh_taper_alternatives.json')]
    if not all(os.path.exists(q) for q in srcs):
        skip('%s length note against its computations' % name, 1, 'provenance/rc6_scripts outputs missing')
        continue
    lb = SP.load(os.path.join(ROOT, 'presets', name + '.yaml'))['paper']['length_bound']
    jbf = json.load(open(srcs[0]))
    jb = jbf['bounds']
    sch = {round(q['F_target'], 4): q for q in json.load(open(srcs[1]))[name]['plans']}
    tap = {round(q['F'], 4): q for q in json.load(open(srcs[2]))['readings']}
    bad = []
    for key, val in lb['min_length_mm_constant_D_head'].items():
        got = [b['length_mm_at_pitch'] for b in jb if abs(b['F'] - float(key[2:])) < 1e-4 and 'constant D' in b['counted']]
        if len(got) != 1 or abs(got[0] - val) > 0.005:
            bad.append('constant-D bound %s %s vs %s' % (key, val, got))
    for key, val in lb['min_length_mm_this_head'].items():
        got = [b['length_mm_at_pitch'] for b in jb if abs(b['F'] - float(key[2:])) < 1e-4 and 'own head' in b['counted']]
        if len(got) != 1 or abs(got[0] - val) > 0.005:
            bad.append('own-head bound %s %s vs %s' % (key, val, got))
    for key, val in lb['this_schedule_head_plus_rows_mm'].items():
        q = sch.get(round(float(key[2:]), 4))
        if q is None or abs(q['head_plus_rows_mm'] - val) > 0.005:
            bad.append('schedule %s %s vs %s' % (key, val, q and q['head_plus_rows_mm']))
    for key, val in lb['linear_taper_rows_mm'].items():
        q = tap.get(round(float(key[2:]), 4))
        if q is None or abs(q['length_mm_linear_75_to_52'] - val) > 0.005:
            bad.append('linear taper %s %s vs %s' % (key, val, q and q['length_mm_linear_75_to_52']))
    q = tap.get(0.505)
    if q is None or abs(q['constant_f_needed_in_printed_length_over_f_max'] - lb['mean_f_needed_over_f_max_at_F_0.505']) > 5e-4:
        bad.append('mean extraction needed at F 0.505')
    sens = json.load(open(srcs[2])).get('sensitivity', [])
    shares = lb.get('largest_share_tapering_at_F_0.505', {})
    sp_, _ = SP.to_engine(SP.load(os.path.join(ROOT, 'presets', name + '.yaml')))
    pitch_ = sp_.G + sp_.g_p
    if not (abs(lb.get('pitch_um', -1) - pitch_) < 1e-9 and abs(jbf.get('pitch_um', -1) - pitch_) < 1e-9 and
            abs(lb.get('pitch_tolerance_um', -1) - dn['post_pitch_um']['tolerance']) < 1e-9):
        bad.append('pitch %s +- %s (note) / %s (bound) vs the preset\'s G + g_p %s and the measured tolerance %s' % (
            lb.get('pitch_um'), lb.get('pitch_tolerance_um'), jbf.get('pitch_um'), pitch_, dn['post_pitch_um']['tolerance']))
    for key, val in shares.items():
        hits = [q for q in sens if abs(q['F'] - 0.505) < 1e-9 and abs(q['pitch_um'] - float(key.split('_')[1])) < 1e-6]
        r_ = hits[0] if len(hits) == 1 else None
        got = None if r_ is None else (r_['largest_share_tapering_with_%s_head_credit' % name.split('_')[-1]]
                                       if key.endswith('head_credit') else r_['largest_share_tapering'])
        if r_ is None or got != val:
            bad.append('share tapering %s %s vs %s' % (key, val, got))
    if not shares:
        bad.append('largest_share_tapering_at_F_0.505 missing')
    check('%s: paper.length_bound equals the computations it quotes (%d values)' % (
        name, sum(len(lb[k]) for k in ('min_length_mm_constant_D_head', 'min_length_mm_this_head',
                                        'this_schedule_head_plus_rows_mm', 'linear_taper_rows_mm',
                                        'largest_share_tapering_at_F_0.505')) + 2),
          not bad, '; '.join(bad[:3]))
n_not = sum(n for _, n, _, _ in not_run)
n_unrequested = sum(n for _, n, _, req in not_run if not req)
if fails:
    print('FAILED: %s (%d checks run, %d not run)' % (fails, n_run[0], n_not))
    sys.exit(1)
if n_unrequested:
    print('INCOMPLETE: %d checks ok, %d not run (%s); not a pass.  Build the presets (runs/ with its pickles) for the '
          'full suite, or run with --source-anchors-only for the checks that need no build'
          % (n_run[0], n_not, '; '.join('%s: %s' % (w, y) for w, _, y, _ in not_run)))
    sys.exit(4)
print('ALL OK (%d checks%s)' % (n_run[0], '; source anchors only: %d drawn-mask and record checks not run, by request'
                                % n_not if n_not else ''))
sys.exit(0)
