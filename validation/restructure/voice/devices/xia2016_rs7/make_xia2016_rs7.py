"""make_xia2016_rs7.py: device 4 of 4 of the opencif 1.0.0 device kits, the RS-7 leukocyte-reduction device of Xia et
al. 2016 (Sci. Rep. 6:35943), reproduced with opencif and aligned with Fig. 3.

    python3 make_xia2016_rs7.py            all steps: check, build, trace, measure, render
    python3 make_xia2016_rs7.py render     one step (every run re-checks the code and the inputs first)

The kit builds two photomasks with the same code, the release's opencif under tree/:
  aligned    inputs/xia2016_rs7_aligned.yaml, rc6's traced variant xia2016_rs7_depicted_fig3 with the parallel-plate
             taper that Fig. 3's panels d and e read, and with both targets fitted to Fig. 3 at that taper (F 0.94072
             from panel g's last side channel, 624.6 um; D 6.415908 um, the cutoff at which the rows fill Fig. 3a's
             legs at that F). The taper is `centre_channel.shear_formula: parallel_plates`, a reproduce-mode departure
             from spec v5 S16-b. It is a reading of the figures and has not been shown to be Xia's method, because
             Xia's text does not state the shear coefficient. It is the package's device, and its figures are
             out/fig_*.
  printed    inputs/xia2016_rs7_depicted.yaml, the printed targets (f_ref, ~15 : 1) and rc6's exact constant-shear
             taper, which is the reference (out/printed/).

The steps run in this order, and every step stops with exit code 2 if a check fails:
  check    the package's opencif against the release (its code stamp against tree/validation/RELEASE.json) and every
           file of the package against MANIFEST.json, the trace code among them; the inputs' sha256 against
           reference/EXPECT.json, which holds RC8's records under the release's code stamp; and each preset against
           its explicit table (INPUTS) in both directions
  build    both photomasks at once, each in a Python process of its own and each byte-identical to the release's
           photomask of the same preset (which is RC8's rebuild), with its exit code. The step then runs the frozen
           physics' own checks of the drawn device (network F = target, no undeclared row outside the cutoff band, one
           fluid component, nothing lost), with the departures that the build declares for this device.
  trace    for each photomask, the trace scripts, run as RC8 ran them (provenance/rc7_traces/code/run_all_rc7.sh).
           They register the whole device on Fig. 3a (xia_fig3a_device.py), panel g (run_xia.py), and panels c, d,
           e, and f (run_xia_more.py, which places c and d by their own features and e and f at Fig. 3a's box
           centers, so that their V12 is partly by construction), and then Fig. 3a again with the traced panels
           (V12). Every record must equal RC8's (reference/), and the trace scripts' output is kept in
           out/trace_log_<photomask>.txt.
  measure  the alignment of each photomask with the figures, in out/alignment.json. It holds Fig. 3a's group shifts
           (V11) and boxes (V12, each logged); each panel's metric at the trace's own window (at most half the gap,
           rc7's window rule), its null ratios, and its posts' flanks along the flow; the physics of each drawn device
           (the share of F within 1% of D, the declared departures, the reverse flow, and the plateau); and the
           panel-f open item (register INV-4), which is panels f's and g's side-channel walls against the mask's,
           measured here by xia_side_walls.py.
  render   the aligned device's figures, out/fig_zoom_overlay, fig_zoom_micrographs (panels c, d, e, f, g),
           fig_full_device, fig_full_device_stacked (the micrographs placed by their registrations), and
           fig_drawing_overlay (Fig. 3a, the paper's whole-device drawing, with the mask placed by its top-anchored
           pose); the printed reference's zoom overlay (out/printed/fig_zoom_overlay); and out/registration.json,
           out/run_log.txt (every run's lines), and out/SHA256SUMS

The figures add no text, so that the micrographs keep the publisher's own labels. Every number is in
out/registration.json and out/alignment.json, and README.md is written from them and from the build records."""
import json
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import devicekit as K

K.SKIP_KEYS |= {'dxf', 'record'}           # Fig. 3a's record names its input files (paths, skipped as strings)
RASTERS = os.path.join(HERE, 'inputs')
OUT = os.path.join(HERE, 'out')
REF = os.path.join(HERE, 'reference')
IMAGE = os.path.join(RASTERS, 'xia5-000.png')
with open(os.path.join(REF, 'EXPECT.json')) as _fh:
    EXPECT = json.load(_fh)
DEV = {
    'aligned': dict(name='xia2016_rs7_aligned', tag='pp', out=OUT, label='the mask xia2016_rs7_aligned'),
    'printed': dict(name='xia2016_rs7_depicted', tag='rc7', out=os.path.join(OUT, 'printed'),
                    label='the mask xia2016_rs7_depicted'),
}
for d in DEV.values():
    d['preset'] = os.path.join(RASTERS, d['name'] + '.yaml')
    d['photomask'] = os.path.join(d['out'], 'photomask')
    d['dxf'] = os.path.join(d['photomask'], d['name'] + '.dxf')
    d['record'] = os.path.join(d['photomask'], d['name'] + '.json')
    d['traces'] = os.path.join(d['out'], 'registrations')
    e = EXPECT['builds'][d['name']]
    d['dxf_sha256'], d['exit'], d['declared'] = e['dxf_sha256'], e['exit_code'], tuple(e['declared'])
    with open(d['preset']) as _fh:
        _p = K.yaml_load(_fh)['process']
    d['post_max'] = (float(_p['post_length_um']) + 1.0, float(_p['post_width_um']) + 1.0)   # the trace's post filter
TRACE_FILES = ['xia_fig3a_device_%s.json', 'xia_results_%s.json'] + ['xia_results_fig3%s_%%s.json' % p for p in 'cdef']

FILES = {os.path.join(HERE, k): v for k, v in EXPECT['inputs_sha256'].items()}

# The generator inputs, written out. The codes are P, printed in the paper; M, measured on its figures; D, derived from
# printed values through the frozen model (or fitted to a measurement, which the entry says); A, assumed where the
# paper is silent; and G, the generator's own choice (docs/PRESET_PROVENANCE.md).
_COMMON = [
    ('process.depth_um', 140, 'A', 'printed 135-150 um'),
    ('process.gap_um', 16, 'P', 'Fig. 1b(ii): G = 16 um'),
    ('process.post_length_um', 35.8, 'M', 'g_p (along the flow): the row pitch 51.78 um on Figs. 3c-g at 1.798 um/px minus G'),
    ('process.post_width_um', 70.9, 'M', 'l_g (across the flow): the post columns\' separation minus the printed centre (Figs. 3c, 3g)'),
    ('process.min_feature_um', 16, 'P', 'Xia line 99: the minimum feature "is now equivalent to its gap size, G"'),
    ('process.min_radius_um', 4.5, 'G', 'the corner limit'),
    ('centre_channel.taper', 'constant_shear', 'A', 'Xia prints the principle (a relatively constant maximum shear) and the ends, not the profile: opencif\'s constant-shear form (LIMITATIONS 16, 24 c)'),
    ('centre_channel.width_in_um', 300, 'P', 'Fig. 1b(i): w_c(i = 0) = 300 um'),
    ('centre_channel.width_min_um', 150, 'P', 'Fig. 1d(i): w_c decreases until it reaches 150 um'),
    ('layout.style', 'paperclip', 'M', 'Fig. 3a: nested legs'),
    ('layout.legs', 7, 'P', 'Fig. 3 caption: "the seventh leg of the device"; Fig. 3a'),
    ('meander.enabled', True, 'P', 'Fig. 1b: the centre channel flanked by serpentine side segments'),
    ('meander.style', 'hairpins', 'M', 'Fig. 1b, Fig. 3b: serpentine (hairpin) segments'),
    ('meander.channel_um', 22.8, 'P', 'Fig. 1b(ii): GS = 22.8 um'),
    ('meander.wall_um', 35.8, 'M', 'Fig. 3b\'s unit pitch and Fig. 3c\'s block pitch, 117.1 um = 2 (GS + wall)'),
    ('meander.head_width_um', 1840, 'M', 'the head envelope on Fig. 3 at 1.798 um/px'),
    ('meander.seal', 'allow', 'G', 'the generator may seal gap 1 with the head'),
    ('meander.transition', 'paper', 'M', 'Fig. 1c: the papers\' rectilinear transition'),
    # `meander.rectilinear_rows` (9, Fig. 1c) is no longer in these presets, because their depicted head
    # (`reproduction.depicted_head`) draws its own 76 rows, and rc7 refuses a spec value that a reproduction field
    # replaces (register BUILD-2).
    ('meander.seam_side_width_um', 16, 'P', 'Fig. 1c(ii): w_s(i = 1) equals the gap'),
    ('meander.seam_at', 'post', 'G', 'where the head meets the post rows'),
    ('inlet.style', 'plain', 'M', 'Fig. 3a'),
    ('inlet.port_radius_um', 862, 'M', 'Fig. 3a at its 5 mm bar'),
    ('inlet.stem_half_width_um', 166, 'M', 'Fig. 3a'),
    ('inlet.port_gap_um', 591, 'M', 'Fig. 3a'),
    ('outlet.style', 'ring', 'M', 'Fig. 3a: a ring outlet'),
    ('outlet.port_radius_um', 1718, 'M', 'Fig. 3a'),
    ('outlet.centre_port_radius_um', 1723, 'M', 'Fig. 3a'),
    ('outlet.tubing.centre_id_mm', 0.25, 'P', 'p. 10'),
    ('outlet.tubing.side_id_mm', 0.86, 'P', 'p. 10'),
    ('outlet.tubing.length_m', 0.15, 'A', 'only the I.D.s are printed'),
    ('outlet.ring.radius_um', 4195, 'M', 'Fig. 3a'),
    ('outlet.ring.centre_offset_um', 1191, 'M', 'Fig. 3a'),
    ('outlet.ring.filtrate_port_offset_um', 4781, 'M', 'Fig. 3a'),
    ('outlet.ring.turn_radius_um', 250, 'G', 'the turn from the side channel onto the ring'),
    ('gates.mode', 'reproduce', 'G', 'draw what the paper depicts and report the envelope'),
    ('reproduction.margin', 20.0, 'G', 'the head\'s fitting margin'),
    ('reproduction.ws_rect', 22.8, 'M', 'Fig. 1c: the side channel beside the transition\'s blocks, GS wide'),
    ('reproduction.leg_Y', [0.0, 1332.0, -1745.0, 2229.0, -2650.0, 3267.0, -3975.0], 'M', 'Fig. 3a: the seven leg axes'),
    ('reproduction.pc_inlet_offset', 2553.0, 'M', 'Fig. 3a: where leg 1 starts below the inlet'),
    ('reproduction.pc_out_ext', 3441.0, 'M', 'Fig. 3a: the outlet stretch'),
    ('reproduction.loop_pull_in', 2.0, 'A', 'Xia: "subtract 1-5 um from the predicted side channel widths"'),
    ('reproduction.loop_pull_in_leg', True, 'P', 'Xia: "of the loops themselves, as well as the downstream side channels of the next device leg"'),
    ('reproduction.loop_ramp', 40.0, 'G', 'the pull-in\'s ramp length'),
    ('reproduction.depicted_head', 'xia', 'M', 'the head as depicted on Figs. 3b and 3a (10 units measured, the rest extrapolated: A)'),
    ('operation.flow_uL_min', 800, 'A', 'inside the printed range ("~1 mL/min", ">0.75 mL/min")'),
]
ALIGNED_D = 6.415908497774869   # D (fitted): provenance/rc7_scripts/xia_fig3_fit_pp.json
ALIGNED_F = 0.9407199784857051  # D (fitted): ibid.
INPUTS = {
    'printed': _COMMON + [
        ('targets.cutoff_um', 6.8644008043574924, 'D', 'printed f_ref 5.8e-4 at w_c(ref) 100 um (corrected eqn 1) through kappa_g'),
        ('targets.filtrate_fraction', 0.938, 'D', 'printed ~15:1 filtrate:retentate (Fig. 3 caption)')],
    'aligned': _COMMON + [
        ('targets.cutoff_um', ALIGNED_D, 'D', 'fitted: the cutoff at which the rows fill Fig. 3a\'s legs at this F, with the parallel-plate taper (provenance/rc7_scripts/xia_fig3_fit_pp.py; x 0.9347 of the printed f_ref\'s)'),
        ('targets.filtrate_fraction', ALIGNED_F, 'D', 'fitted: the F at which the last side channel is Fig. 3g\'s 624.6 um (15.9:1; printed ~15:1)'),
        ('centre_channel.shear_formula', 'parallel_plates', 'A', 'rc7 (XIA-1): the taper Fig. 3\'s panels d and e read, chosen with them in view; a reproduce-mode departure from spec v5 S16-b (Xia\'s text does not state the coefficient)')],
}
PANELS = ['c', 'd', 'e', 'f', 'g']
# The panels that are recorded rather than claimed (the trace anchors' role, provenance/rc7_traces/trace_anchors.json).
# As in rc7, they are panel e on the printed (depicted) mask and g on the aligned one, because their two registrations
# disagree (V1), and f on both, because it is read at its own scale (register INV-4). RC8 added d on both masks and g on
# the printed one, because the metric covers under half of their eligible samples and their censored null ratios are
# near 1 (D5 of the adjudication of the rc7 reviews; A7). The claims are c on both masks and e on the aligned one.
RECORDED = {('aligned', 'd'), ('printed', 'd'), ('printed', 'g'), ('printed', 'e'), ('aligned', 'g'), ('aligned', 'f'),
            ('printed', 'f')}
KEYS = {'c': 'xia_results_fig3c_%s.json', 'd': 'xia_results_fig3d_%s.json', 'e': 'xia_results_fig3e_%s.json',
        'f': 'xia_results_fig3f_%s.json', 'g': 'xia_results_%s.json'}
STACK_ROWS = 3     # the whole device in three rows (its footprint is the build record's drawn_parameters.footprint_mm)
PLACED_BY = dict(c='the first pill below a non-periodic feature (the leg\'s first post)',
                 d='the first pill below a non-periodic feature on each leg (the leg\'s first post)',
                 e='Fig. 3a\'s box centre (its V12 partly by construction)', f='Fig. 3a\'s box centre (its V12 partly by construction)',
                 g='leg 7\'s last posts, indexed from the last one (the outlet)')


def step_check():
    K.say('== check')
    K.check_code()
    K.say('reference/EXPECT.json (the frozen tree\'s records): sha256 %s' % K.sha256(os.path.join(REF, 'EXPECT.json'))[:16])
    if EXPECT['code_stamp'] != K.STAMP:
        K.fail('reference/EXPECT.json is for code %s, the package for %s' % (EXPECT['code_stamp'], K.STAMP))
    K.check_inputs(FILES)
    for key in ('aligned', 'printed'):
        K.say('generator input (%s):' % key)
        K.check_preset(DEV[key]['preset'], INPUTS[key])


def _build_one(key):
    """Build one photomask in this process, which step_build starts as a Python process of its own."""
    d = DEV[key]
    K.say('build (%s), in a Python process of its own:' % key)
    K.build(d['preset'], d['photomask'], d['dxf_sha256'], d['exit'])


def step_build():
    """Build each photomask in a Python process of its own, the two at once, as RC8's rebuild built them. Because
    pipeline.run clears the memo caches first (BUILD-1), one process would do as well, and the separate processes also
    run the two builds at once."""
    K.say('== build')
    jobs = {}
    for key in ('aligned', 'printed'):
        os.makedirs(DEV[key]['photomask'], exist_ok=True)
        log = os.path.join(OUT, '.build_%s.txt' % key)
        jobs[key] = (log, subprocess.Popen([sys.executable, '-I', '-B', os.path.abspath(__file__), '_build_one', key],
                                           stdout=open(log, 'w'), stderr=subprocess.STDOUT))
    for key, (log, p) in jobs.items():
        p.wait()
        text = open(log).read()
        os.remove(log)
        sys.stdout.write(text)
        sys.stdout.flush()
        K.LOG.extend(l for l in text.splitlines() if re.match(r'\d{4}-\d\d-\d\d \d\d:\d\d:\d\d  ', l))
        if p.returncode != 0:
            K.fail('the %s build stopped (exit %d; its lines are above)' % (key, p.returncode))
    for key in ('aligned', 'printed'):
        d = DEV[key]
        K.say('physics (%s):' % key)
        K.check_build_physics(json.load(open(d['record'])), d['declared'])


def _need(key):
    """Return the photomask's entry, and stop unless out/ holds its photomask with the reference sha256."""
    d = DEV[key]
    if not os.path.exists(d['dxf']) or K.sha256(d['dxf']) != d['dxf_sha256']:
        K.fail('out/ holds no %s photomask matching its reference: run the build step' % key)
    return d


def _run(args, log):
    """Run one of the trace scripts as RC8 ran it, with python -I so that only the package's own code is imported,
    append its output to log, and return its run time in seconds."""
    t0 = time.time()
    p = subprocess.run([sys.executable, '-I'] + args, capture_output=True, text=True)
    open(log, 'a').write('$ %s\n%s\n%s\n' % (' '.join(args), p.stdout, p.stderr))
    if p.returncode != 0:
        K.fail('%s exited %d (its output: %s)' % (os.path.basename(args[0]), p.returncode, log))
    return time.time() - t0


def step_trace():
    K.say('== trace')
    code = K.TRACE
    for key in ('aligned', 'printed'):
        d = _need(key)
        tag = d['tag']
        work = os.path.join(OUT, 'trace_work')            # the trace scripts' captioned figures (raster crops), which
        os.makedirs(work, exist_ok=True)                  # are removed afterwards, whatever happens
        log = os.path.join(OUT, 'trace_log_%s.txt' % key)  # the trace scripts' own output, kept (package check, 28)
        try:
            f3a = os.path.join(work, 'xia_fig3a_device_%s.json' % tag)
            common = [d['preset'], d['dxf'], d['record'], d['label']]
            t = _run([os.path.join(code, 'xia_fig3a_device.py'), IMAGE, d['dxf'], d['record'], f3a,
                      os.path.join(work, 'xia_fig3a_device_%s.png' % tag)], log)
            K.say('trace (%s): the whole device on Fig. 3a, %.0f s' % (key, t))
            t = _run([os.path.join(code, 'run_xia.py'), K.TREE, RASTERS, work, tag, 'rc6'] + common, log)
            K.say('trace (%s): panel g, %.0f s' % (key, t))
            t = _run([os.path.join(code, 'run_xia_more.py'), K.TREE, RASTERS, work, tag, 'cdef'] + common + ['--fig3a', f3a], log)
            K.say('trace (%s): panels c, d, e, f, %.0f s' % (key, t))
            traced = ['%s:-:%s' % (os.path.join(work, KEYS[p] % tag), p) for p in 'cdefg']
            t = _run([os.path.join(code, 'xia_fig3a_device.py'), IMAGE, d['dxf'], d['record'], f3a,
                      os.path.join(work, 'xia_fig3a_device_%s.png' % tag)] + traced, log)
            K.say('trace (%s): Fig. 3a with the traced panels (V12), %.0f s' % (key, t))
            os.makedirs(d['traces'], exist_ok=True)
            for pat in TRACE_FILES:
                fn = pat % tag
                K.check_trace(json.load(open(os.path.join(work, fn))), os.path.join(REF, fn), os.path.join(d['traces'], fn))
        finally:
            shutil.rmtree(work, ignore_errors=True)


def step_measure():
    K.say('== measure')
    K.use_tree()
    import trace_lib as T
    import gifford_alignment as XA                    # for narrow_window, the same function that the Gifford kits use
    import xia_side_walls as XW                       # the trace code's panel f / g wall reading (register INV-4)
    raster = T.Raster(IMAGE, EXPECT['inputs_sha256']['inputs/xia5-000.png'][:16])
    imf = raster.im.astype(float)
    res = dict(not_blind=('fitted to the figures, so consistency rather than tests: F (to panel g\'s 624.6 um last side '
                          'channel), D (to Fig. 3a\'s leg fill) and the taper (chosen with panels d and e in view); the '
                          '1.798 um/px scale is in part calibrated on panels c\'s and g\'s printed widths; blind: the bars, '
                          'the rotation (the image) and V1 (two independent registrations agree)'))
    for key in ('aligned', 'printed'):
        d = _need(key)
        mask = T.Mask(d['dxf'])
        f3a = json.load(open(os.path.join(d['traces'], 'xia_fig3a_device_%s.json' % d['tag'])))
        r = dict(fig3a_group_shifts_um={g: dict(along=v['shift_along_um'], across=v['shift_across_um'])
                                        for g, v in f3a['V11_group_shifts'].items()},
                 fig3a_boxes={b: dict(ok=v.get('ok'), offset_over_half_size=v.get('offset_over_half_size'))
                              for b, v in f3a['V12_boxes'].items() if 'traced' in v})
        for g, v in r['fig3a_group_shifts_um'].items():      # None along a leg's straight part (no along-flow feature)
            K.say('measure (%s): Fig. 3a, %-24s along %8s um, across %+6.0f um' % (
                key, g, '%+.0f' % v['along'] if v['along'] is not None else '-', v['across']))
        for b, v in sorted(r['fig3a_boxes'].items()):
            K.say('measure (%s): Fig. 3a, panel %s %s its box (%.2f of the box\'s half-size; placed by %s)' % (
                key, b, 'inside' if v['ok'] else 'OUTSIDE', v['offset_over_half_size'], PLACED_BY.get(b, '?')))
        panels = {}
        for p in PANELS:
            rec = json.load(open(os.path.join(d['traces'], KEYS[p] % d['tag'])))
            wpx = float(rec['metric_search_px'])
            nw = XA.narrow_window(T, raster, mask, rec, d['post_max'], lambda e: True, windows=(wpx,))
            w_ = nw['%.1f' % wpx]
            nr = K.null_ratio(rec)
            panels[p] = dict(metric=rec['metric'], narrow_window=nw, measured=rec['measured'], placed_by=PLACED_BY[p],
                             metric_window_um=rec.get('metric_window_um'), window_rule=rec.get('window_rule'),
                             V1=rec['V1_feature_vs_chamfer_px'], null_ratio=nr,
                             below_null_criterion=sorted(k for k, v in nr.items() if v < K.NULL_RATIO_CRITERION),
                             coverage=rec['metric'].get('coverage'), n_eligible=rec['metric'].get('n_eligible'),
                             status='recorded' if (key, p) in RECORDED else 'claimed')
            K.say('measure (%s): panel %s: median / p90 %.2f / %.2f um at the trace\'s +-%.1f px window (%.1f um, at '
                  'most half the gap), coverage %.3f (%s); posts\' flanks along the flow %+.2f um; V1 %s; null ratios %s' % (
                      key, p, rec['metric']['median_abs_um'], rec['metric']['p90_abs_um'], wpx, w_['window_um'],
                      rec['metric'].get('coverage') or float('nan'), panels[p]['status'],
                      w_['post_flanks_along']['median_um'], 'agree' if rec['V1_feature_vs_chamfer_px']['agree'] else 'DISAGREE',
                      ', '.join('%s %.2f%s' % (k, v, ' (below %.1f)' % K.NULL_RATIO_CRITERION if v < K.NULL_RATIO_CRITERION else '')
                                for k, v in sorted(nr.items()))))
            if p in 'fg':
                sw = XW.side_walls(imf, rec, mask)
                panels[p]['side_walls'] = sw
                K.say('measure (%s): panel %s: the side channels\' outer walls %s (image inside the mask\'s, at %.3f um/px)' % (
                    key, p, ', '.join('%s %+.1f um' % (w['side'], w['image_inside_um']) for w in sw['walls']), sw['um_per_px']))
        r['panels'] = panels
        r['physics'] = K.physics_block(json.load(open(d['record'])))
        res[key] = r
    res['open_item'] = ('panel f (a quarter of the way down leg 7): the image\'s side-channel outer walls lie inside the '
                        'mask\'s (panels.f.side_walls, at panel f\'s own pitch scale; about 7 um less at the 1.798 um/px '
                        'of the other panels), on either taper, while panel g\'s wall agrees: register INV-4, open')
    with open(os.path.join(OUT, 'alignment.json'), 'w') as fh:
        json.dump(res, fh, indent=1, default=K.jsonable)
    K.say('measure: -> out/alignment.json')


def _panels(key, T):
    """Return the photomask, its panels, and their registration records, after checking that every panel was
    registered on this mask and claims every edge in its frame."""
    d = _need(key)
    mask = T.Mask(d['dxf'])
    recs = {p: json.load(open(os.path.join(d['traces'], KEYS[p] % d['tag']))) for p in PANELS}
    for p, rec in recs.items():
        if rec['mask_sha256'] != mask.sha256:
            K.fail('panel %s was registered on another mask' % p)
        if rec['other_edges']:
            K.fail('panel %s: the trace leaves %d edges unclaimed; this script draws every edge as claimed' % (p, rec['other_edges']))
    rasters = {K.sha256(IMAGE): K.load_rgb(IMAGE)}
    return mask, K.prepare_panels(PANELS, recs, rasters, mask), recs


def step_render():
    K.say('== render')
    K.use_tree()
    import trace_lib as T
    mask, panels, recs = _panels('aligned', T)
    z1 = K.render_zoom(panels, mask, os.path.join(OUT, 'fig_zoom_overlay'), overlay=True)
    z0 = K.render_zoom(panels, mask, os.path.join(OUT, 'fig_zoom_micrographs'), overlay=False)
    K.say('render: zoomed panels c, d, e, f, g: %d x %d px (each image pixel x %d), with and without the overlay' % (
        z1['width_px'], z1['height_px'], z1['magnify']))
    full = K.render_full(mask, panels, os.path.join(OUT, 'fig_full_device'))
    K.say('render: whole device: PDF page %.1f x %.1f in (%g x the device), PNG %d x %d px at %g um/px' % (
        full['pdf']['page_in'][0], full['pdf']['page_in'][1], full['pdf']['magnification'],
        full['png']['size_px'][0], full['png']['size_px'][1], full['png']['um_per_px']))
    stack = K.render_stack(mask, panels, os.path.join(OUT, 'fig_full_device_stacked'), full, rows=STACK_ROWS)
    K.say('render: whole device in %d stacked rows of %s mm (micrographs: %s): PDF page %.1f x %.1f in, '
          'PNG %d x %d px' % (stack['rows'], '/'.join('%.2f' % v for v in stack['row_length_mm']),
                              ', '.join('%s in row %d' % kv for kv in stack['panel_row'].items()),
                              stack['pdf']['page_in'][0], stack['pdf']['page_in'][1], *stack['png']['size_px']))
    f3a = json.load(open(os.path.join(DEV['aligned']['traces'], 'xia_fig3a_device_%s.json' % DEV['aligned']['tag'])))
    und = dict(f3a['underlay'], micrograph_box=f3a['underlay']['crop'])
    rasters = {K.sha256(IMAGE): K.load_rgb(IMAGE)}
    dpan = K.prepare_panels(['a'], {'a': und}, rasters, mask)
    zd = K.render_zoom(dpan, mask, os.path.join(OUT, 'fig_drawing_overlay'), overlay=True)
    K.say('render: Fig. 3a (the whole-device drawing) with the mask placed by its top-anchored pose: %d x %d px' % (
        zd['width_px'], zd['height_px']))
    mask_p, panels_p, _ = _panels('printed', T)
    zp = K.render_zoom(panels_p, mask_p, os.path.join(DEV['printed']['out'], 'fig_zoom_overlay'), overlay=True)
    K.say('render: the printed reference\'s zoomed panels with their overlay, %d x %d px' % (zp['width_px'], zp['height_px']))
    summary = dict(device=DEV['aligned']['name'], paper='Xia et al. 2016, Sci. Rep. 6:35943, Fig. 3 (panels c-g; panel a)',
                   code_stamp=K.STAMP,
                   photomask=dict(file='photomask/%s.dxf' % DEV['aligned']['name'], sha256=mask.sha256),
                   zoom_overlay=z1, zoom_micrographs=z0, full_device=full, full_device_stacked=stack, drawing_overlay=zd,
                   drawing_pose=dict(underlay=f3a['underlay'], pose_top=f3a['pose_top']),
                   panels=K.registration_summary(panels),
                   reference_printed=dict(photomask=dict(file='printed/photomask/%s.dxf' % DEV['printed']['name'],
                                                         sha256=mask_p.sha256), zoom_overlay=zp,
                                          panels=K.registration_summary(panels_p)),
                   colours=dict(cyan='the mask\'s edges (every panel claims all the edges in its frame)'),
                   coordinates='mask um: x along the device, y across; image px: (column, row) of the raster; '
                               'image px = origin + x ex + y ey (map_chamfer)')
    with open(os.path.join(OUT, 'registration.json'), 'w') as fh:
        json.dump(summary, fh, indent=1)


STEPS = dict(check=step_check, build=step_build, trace=step_trace, measure=step_measure, render=step_render)


def main(argv):
    todo = argv or list(STEPS)
    os.makedirs(OUT, exist_ok=True)
    try:
        if 'check' not in todo:
            step_check()
        for s in todo:
            STEPS[s]()
        K.say('done: %s' % ', '.join(todo))
    except Exception as e:
        K.say('STOPPED: %s: %s' % (type(e).__name__, e))
        raise
    finally:
        with open(os.path.join(OUT, 'run_log.txt'), 'a') as f:      # every run's lines kept (package check, 28)
            f.write('\n'.join(K.LOG) + '\n')
        if 'render' in todo:
            K.write_sums(OUT)


if __name__ == '__main__':
    if sys.argv[1:2] == ['_build_one']:          # the build step's own call, for one photomask and nothing else
        _build_one(sys.argv[2])
    else:
        main(sys.argv[1:])
