"""make_gifford2014_fig6.py: device 3 of 4 of the opencif 1.0.0 device kits, the Fig. 6 platelet device of Gifford et
al. 2014 (Lab Chip 14, 4496-4505), reproduced with opencif and aligned with the Fig. 6 micrograph.

    python3 make_gifford2014_fig6.py            all steps: check, build, trace, measure, render
    python3 make_gifford2014_fig6.py render     one step (every run re-checks the code and the inputs first)

The kit builds two photomasks with the same code, the release's opencif under tree/:
  aligned    inputs/gifford2014_fig6_aligned.yaml, Gifford's printed design with the text's final side width (144 um,
             where the caption says 140), his start rule at a 40 um first side channel, and the micrograph's 17.1 um
             gap. The start width is fitted on leg 1 (D, fitted), and legs 3 and 5 test it only in part, because the
             text's 144 um was preferred with them in view; leg 7's end is set by F. It is the package's device, and
             its figures are out/fig_*.
  design     inputs/gifford2014_fig6_design.yaml, the design-level preset (sealed start, the caption's 140 um, the
             20.5 um gap), which is the reference (out/design/).

The steps run in this order, and every step stops with exit code 2 if a check fails:
  check    the package's opencif against the release (its code stamp against tree/validation/RELEASE.json) and every
           file of the package against MANIFEST.json, among them the trace code and the two measurement files that
           the Fig. 6 trace reads; the inputs' sha256 against reference/EXPECT.json, which holds RC8's records under
           the release's code stamp; and each preset against its table in both directions
  build    both photomasks, each byte-identical to the release's photomask of the same preset (which is RC8's
           rebuild) and with its exit code (that of RC8's record). The step runs the frozen physics' own checks of the
           drawn device (network F = target, no undeclared row outside the cutoff band, one fluid component, nothing
           lost) and writes its physics block (the plateau, the cutoff band below it, the declared departures, and
           the model band's junction terms). The aligned device declares its start rows ('floor start').
  trace    the micrograph registered on each photomask by the trace code, as RC8 traced it, with each of legs 7, 5,
           3, and 1 on its own map and the four legs with their turns on one map. Every value must equal RC8's
           records (reference/, from provenance/rc7_traces/results/gifford/).
  measure  (1) the start width that the micrograph implies. Leg 1's side channels on the image are compared with the
           sealed reference's drawn width there, through the row-fraction law (gifford_alignment.py), at the text's
           final width. The width is fitted on leg 1, legs 3 and 5 test it in part, and the aligned preset's start
           must be that value to the micrometer. (2) the alignment of the aligned photomask, with every leg-end side
           channel within one image pixel, every leg's median below 1 um, the posts' flanks along the flow within
           0.5 um, the turns' walls within 2 um, and leg 7's drawn end beyond the image's (register GIF6-4)
           -> out/alignment.json. These are regression limits, set with the results in view.
  render   the aligned device's figures (out/fig_zoom_overlay, fig_zoom_micrographs, fig_full_device,
           fig_full_device_stacked; PDF + PNG each), the reference's zoom overlay (out/design/fig_zoom_overlay),
           out/registration.json, out/run_log.txt (every run's lines), and out/SHA256SUMS

The figures add no text, so that the micrograph keeps the publisher's own labels. Every number is in
out/registration.json and out/alignment.json, and README.md is written from them and from the build records."""
import json
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import devicekit as K

RASTERS = os.path.join(HERE, 'inputs')
OUT = os.path.join(HERE, 'out')
REF = os.path.join(HERE, 'reference')
IMAGE = os.path.join(RASTERS, 'gif7-000.png')
with open(os.path.join(REF, 'EXPECT.json')) as _fh:
    EXPECT = json.load(_fh)
DEV = {
    'aligned': dict(name='gifford2014_fig6_aligned', out=OUT, tag='aligned', label='the mask gifford2014_fig6_aligned'),
    'design': dict(name='gifford2014_fig6_design', out=os.path.join(OUT, 'design'), tag='design',
                   label='the mask gifford2014_fig6_design'),
}
for d in DEV.values():
    e = EXPECT['builds'][d['name']]
    d.update(exit=e['exit_code'], declared=tuple(e['declared']), dxf_sha256=e['dxf_sha256'])
    d['ref_per_leg'] = os.path.join(REF, 'gifford6_results_%s.json' % d['tag'])
    d['ref_one_map'] = os.path.join(REF, 'gifford6_results_%s_onemap.json' % d['tag'])
    d['preset'] = os.path.join(RASTERS, d['name'] + '.yaml')
    d['photomask'] = os.path.join(d['out'], 'photomask')
    d['dxf'] = os.path.join(d['photomask'], d['name'] + '.dxf')
    d['record'] = os.path.join(d['photomask'], d['name'] + '.json')
    d['reg_per_leg'] = os.path.join(d['out'], 'registration_full_per_leg.json')
    d['reg_one_map'] = os.path.join(d['out'], 'registration_full_one_map.json')
FILES = {os.path.join(HERE, k): v for k, v in EXPECT['inputs_sha256'].items()}

# The generator inputs, written out. The codes are P, printed in the paper; M, measured on its figures; D, derived from
# printed values through the frozen model (or fitted to a measurement, which the entry says); A, assumed where the
# paper is silent; and G, the generator's own choice.
_COMMON = [
    ('targets.cutoff_um', 3.6236, 'D', 'printed f_gap 1.04e-4 at w_c 125 um (Fig. 6 caption), through kappa_g (App. B)'),
    ('process.depth_um', 150, 'P', 'Fig. 6 caption: "The depth of the device was 150 um +-5 um"'),
    ('process.post_width_um', 70.8, 'M', 'l_g (across the flow), on the micrograph at the printed 125 um centre'),
    ('process.min_feature_um', 4.5, 'G', 'the design preset\'s floor (chosen so its sealed start\'s 4.69 um first channel draws)'),
    ('process.min_radius_um', 1.8, 'G', 'ibid.'),
    ('process.fillet_um', 2.0, 'G', 'fillet over the inlet and outlet zones'),
    ('centre_channel.taper', 'constant', 'P', 'Gifford\'s design: a constant central channel'),
    ('centre_channel.width_in_um', 125, 'P', 'Fig. 6 caption: "the width of the central channel was 125 um"'),
    ('layout.style', 'paperclip', 'M', 'Fig. 6 schematic (not drawn to scale, p. 4502): nested legs; the four arrows'),
    ('layout.legs', 7, 'M', 'Fig. 6 schematic: seven legs'),
    ('meander.enabled', False, 'P', 'Gifford\'s devices have no meander'),
    ('inlet.style', 'plain', 'M', 'Fig. 6 schematic: the inlet inside the innermost loop'),
    ('inlet.port_radius_um', 300, 'A', 'the schematic is not to scale'),
    ('outlet.style', 'ring', 'M', 'Fig. 6 schematic: a ring-like outlet (a small centre port inside a ring, one larger port)'),
    ('outlet.port_radius_um', 750, 'A', 'the schematic is not to scale'),
    ('outlet.centre_port_radius_um', 300, 'A', 'ibid.'),
    ('operation.flow_uL_min', 55, 'P', 'p. 4502: "a total throughput of 55 uL min-1"'),
    ('gates.mode', 'reproduce', 'G', 'draw what the paper depicts and report the envelope (below the kappa plateau)'),
    ('reproduction.pc_last_short', 257, 'M', 'the micrograph: leg 7\'s array ends 5 rows (257 um) before legs 5, 3, 1 (reproduction-only)'),
    ('reproduction.leg_Y', [0, 3491.4, -587.8, 4079.3, -1175.4, 4666.8, -1764.3], 'M',
     'legs 1, 3, 5, 7 on the micrograph\'s post rows (588 um apart); legs 2, 4, 6 where semicircular turns about the '
     'turn walls\' common centre put them (the centre measured; the semicircle ASSUMED: 8-23 deg of each turn in frame)'),
]
INPUTS_DESIGN = _COMMON + [
    ('targets.filtrate_fraction', 0.71768, 'D', 'outlet identity at the caption\'s final side width 140 um'),
    ('process.gap_um', 20.5, 'A', 'as Fig. 4 (rc6); this micrograph reads 17.1 um'),
    ('process.post_length_um', 30.9, 'M', 'g_p (along the flow) = the micrograph\'s row pitch 51.38 um minus G 20.5'),
]
INPUTS_ALIGNED = _COMMON + [
    ('targets.filtrate_fraction', 0.72923, 'D', 'outlet identity at the text\'s final side width 144 um (p. 4502; the caption says 140; the micrograph reads 145.5)'),
    ('process.gap_um', 17.1, 'M', 'the micrograph: G_read 17.08 um (rc6 gifford_fig6_measurement.json; Fig. 4 reads 17.2)'),
    ('process.post_length_um', 34.3, 'M', 'g_p (along the flow): the micrograph\'s posts 34.3 um (row pitch 51.38 um minus G 17.1)'),
    ('start.style', 'floor', 'P', 'Gifford\'s start rule (p. 4498), a nonzero first side width: App. A\'s "other start mode" (A:1118-1119); the leg ends imply it (GIF6-1)'),
    ('start.side_width_um', 40.0, 'D', 'fitted: the start that puts leg 1 at the micrograph\'s 58.9 um (the measure step: 40.1 um); legs 3 and 5 test it in part'),
]
PANELS = ['four_legs']    # the registration that the figures use, on one map (legs 7, 5, 3, 1 and the turns)
LEGS = ['leg7', 'leg5', 'leg3', 'leg1']
STACK_ROWS = 3
W_END_TEXT = 144.0        # the text's final side width (p. 4502)


def keyed(recs):
    return {('four_legs' if 'four_legs' in v['name'] else v['name'].split('_')[2]): v for v in recs}


def split_mask(mask, rec, P6, one_map, preset):
    """Rebuild the trace's claim from its record with the trace code's own SplitMask, as rc7's panels_gifford6.run
    builds it, and return (SplitMask, claim(edge) -> bool, post_max). On one map, every edge is claimed except leg 7's
    beyond LEG7_MARGIN_ROWS rows before the image's leg-7 end, which lies LEG7_IMAGE_SHORT_ROWS rows before the mean
    end of legs 5, 3, and 1 (rc7 package check, finding 1). Per leg, each of legs 7, 5, 3, and 1 is claimed to
    CLAIM_PAST_WALL_UM into its divider wall."""
    import yaml
    sp = yaml.safe_load(open(preset))['process']
    g_p, G = sp['post_length_um'], sp['gap_um']
    pitch = g_p + G
    legY = [float(v) for v in rec['mask_record']['leg_Y']]
    x_last = {int(j): float(v) for j, v in rec['measured']['mask_last_post_x_um'].items()}
    half = abs(0.5 * float(K.np.min(K.np.diff(sorted(legY[j - 1] for j in (7, 5, 3, 1))))))
    if one_map:
        x_img7 = sum(x_last[j] for j in (1, 3, 5)) / 3.0 - P6.LEG7_IMAGE_SHORT_ROWS * pitch
        bands = [(legY[6] - half, legY[6] + half, x_img7 - P6.LEG7_MARGIN_ROWS * pitch)]
    else:
        bands = [(legY[j - 1] - half, legY[j - 1] + half, x_last[j] + 0.5 * g_p + G + P6.CLAIM_PAST_WALL_UM)
                 for j in (7, 5, 3, 1)]
    return (P6.SplitMask(mask, bands, outside_claim=one_map), (lambda e: not getattr(e, 'beyond', False)),
            (g_p + 1.0, sp['post_width_um'] + 1.0))


def check_claim(smask, claim, post_max, rec, key):
    """Check that the claim splits the traced region's edges as the trace split them, and return the two counts."""
    es = smask.edges(*rec['region_um'], post_max=post_max)
    n = sum(1 for e in es if claim(e))
    if (n, len(es) - n) != (rec['claim_edges'], rec['other_edges']):
        K.fail('%s: the claim colouring (%d / %d edges) is not the trace\'s (%d / %d)' % (
            key, n, len(es) - n, rec['claim_edges'], rec['other_edges']))
    return n, len(es) - n


def step_check():
    K.say('== check')
    K.check_code()
    K.say('reference/EXPECT.json (the frozen tree\'s records): sha256 %s' % K.sha256(os.path.join(REF, 'EXPECT.json'))[:16])
    if EXPECT['code_stamp'] != K.STAMP:
        K.fail('reference/EXPECT.json is for code %s, the package for %s' % (EXPECT['code_stamp'], K.STAMP))
    K.check_inputs(FILES)
    K.say('generator input (aligned, the package\'s device):')
    K.check_preset(DEV['aligned']['preset'], INPUTS_ALIGNED)
    K.say('generator input (design, the reference):')
    K.check_preset(DEV['design']['preset'], INPUTS_DESIGN)


def step_build():
    K.say('== build')
    for key in ('design', 'aligned'):
        d = DEV[key]
        K.say('build (%s):' % key)
        _, _, rec = K.build(d['preset'], d['photomask'], d['dxf_sha256'], d['exit'])
        K.check_build_physics(rec, d['declared'])


def _need_dxf(key):
    """Return the photomask's entry, and stop unless out/ holds its photomask with the reference sha256."""
    d = DEV[key]
    if not os.path.exists(d['dxf']) or K.sha256(d['dxf']) != d['dxf_sha256']:
        K.fail('out/ holds no %s photomask matching its reference: run the build step' % key)
    return d


def step_trace():
    K.say('== trace')
    K.use_tree()
    import panels_gifford6 as P6
    for key in ('design', 'aligned'):
        d = _need_dxf(key)
        tag = d['tag']
        work = os.path.join(OUT, 'trace_work')            # the trace code's captioned figures (raster crops), which are
        os.makedirs(work, exist_ok=True)                  # removed afterwards, whatever happens
        try:
            for one, ref, reg, what in ((False, d['ref_per_leg'], d['reg_per_leg'], 'legs 7, 5, 3, 1 each on its own map'),
                                        (True, d['ref_one_map'], d['reg_one_map'], 'the four legs and their turns on one map')):
                t0 = time.time()
                recs = P6.run(K.TREE, RASTERS, work, dxf=d['dxf'], record=d['record'], preset=d['preset'], tag=tag,
                              label=d['label'], scan=True, one_map=one, turns=one)
                K.say('trace (%s): Fig. 6, %s, %.0f s' % (key, what, time.time() - t0))
                K.check_trace(keyed(recs), ref, reg)
        finally:
            shutil.rmtree(work, ignore_errors=True)


def step_measure():
    K.say('== measure')
    K.use_tree()
    import gifford_alignment as GA
    res = {}
    pa = K.yaml_load(open(DEV['aligned']['preset']))
    WC, H = float(pa['centre_channel']['width_in_um']), float(pa['process']['depth_um'])

    def legs_table(rec):
        img, drawn = {}, {}
        for r in rec['measured']['side_channels_at_leg_ends']:
            img.setdefault(r['leg'], []).append(r['image_um'])
            drawn.setdefault(r['leg'], []).append(r['mask_um'])
        return img, drawn
    # (1) the start width that the micrograph implies, from the sealed reference's one-map record
    rd = json.load(open(DEV['design']['reg_one_map']))['four_legs']
    img, drawn = legs_table(rd)
    mi = {j: sum(v) / len(v) for j, v in img.items()}
    md = {j: sum(v) / len(v) for j, v in drawn.items()}
    phi_c = {j: GA.phi(md[j], md[7], WC, H) for j in (1, 3, 5)}
    delta, w0, defect = GA.start_from(phi_c[1], mi[1], W_END_TEXT, WC, H)
    pred = {j: GA.predict(phi_c[j], delta, W_END_TEXT, WC, H) for j in (3, 5)}
    w_preset = pa['start']['side_width_um']
    K.say('measure: leg ends on the image %s um; the sealed reference draws %s um (row fractions %s)' % (
        ' / '.join('%.2f' % mi[j] for j in (1, 3, 5, 7)), ' / '.join('%.2f' % md[j] for j in (1, 3, 5, 7)),
        ' / '.join('%.4f' % phi_c[j] for j in (1, 3, 5))))
    K.say('measure: the start leg 1 implies at the text\'s final width %.0f um: %.2f um (the schedule entered at %.2f %% '
          'of its rows; App. A\'s mass defect 2 sigma0 %.3f); it predicts legs 3 and 5 at %.2f and %.2f um against the '
          'image\'s %.2f and %.2f (a partial test: the text\'s 144 um was preferred over the caption\'s 140 with legs 3 '
          'and 5 in view, register GIF6-2); the aligned preset draws %.1f um' % (
              W_END_TEXT, w0, 100 * delta, defect, pred[3], pred[5], mi[3], mi[5], w_preset))
    if abs(round(w0) - w_preset) > 0.01:
        K.fail('the aligned preset\'s start %.1f um is not the measured %.2f um to the micrometre' % (w_preset, w0))
    res['start_width'] = dict(rule='leg 1\'s end: image width against the sealed reference\'s drawn width there, through '
                                   'phi(w) = ln Q(w) / ln Q(w_end) (gifford_alignment.py), at the text\'s final 144 um; '
                                   'legs 3 and 5 test it in part (the text\'s 144 um, which sets leg 7\'s end through F, '
                                   'was preferred over the caption\'s 140 um with legs 3 and 5 in view: register GIF6-2); '
                                   'leg 7 is not a test',
                              image_um=mi, drawn_sealed_um=md, phi_chip=phi_c, delta=delta, start_um=w0,
                              mass_defect_2sigma0=defect, predicted_um=pred, preset_um=w_preset)
    # (2) the alignment of the aligned photomask
    ra1 = json.load(open(DEV['aligned']['reg_one_map']))['four_legs']
    ra4 = json.load(open(DEV['aligned']['reg_per_leg']))
    px_um = float(ra1['measured']['um_per_px'])              # one image pixel, scaled on the printed 125 um center
    checks = []

    def check(what, value, limit, unit='um'):
        ok = abs(value) <= limit
        checks.append(dict(what=what, value=value, limit=limit, ok=ok))
        K.say('measure: %-60s %+8.2f %s (limit %.2f): %s' % (what, value, unit, limit, 'ok' if ok else 'FAILS'))
    for r in ra1['measured']['side_channels_at_leg_ends']:
        check('leg %d side channel (%s side) at its end, image - drawn' % (r['leg'], r['side']), r['image_um'] - r['mask_um'], px_um)
    for key in LEGS:
        check('%s on its own map, median |offset|' % key, ra4[key]['metric']['median_abs_um'], 1.0)
        for c in ('post_along_C', 'post_along_R'):
            v = ra4[key]['V3_by_class'].get(c)
            if v and v['n'] >= 20 and v['mean_signed_um'] is not None:
                check('%s posts\' flanks along the flow (%s), mean offset' % (key, c[-1]), v['mean_signed_um'], 0.5)
    check('one map, median |offset|', ra1['metric']['median_abs_um'], 1.0)
    check('the turns\' walls, median distance to the mask\'s nearest edge', ra1['measured']['turns']['nearest_mask_edge_um']['median'], 2.0)
    le7 = ra1['measured']['leg7_claim_end']
    check('leg 7\'s drawn end beyond the image\'s (rows; GIF6-4)', le7['drawn_end_beyond_image_end_rows'], 4.5, unit='rows')
    res['alignment'] = dict(side_channels_at_leg_ends=ra1['measured']['side_channels_at_leg_ends'], checks=checks,
                            turns=ra1['measured']['turns']['nearest_mask_edge_um'],
                            limits_note=('one image pixel (%.3f um) for widths; 1 um for the medians; 0.5 um for the '
                                         'posts\' flanks; 2 um for the turns (consistency: the axes were read on them); '
                                         '4.5 rows for leg 7\'s drawn end beyond the image\'s (the generator places it '
                                         'in steps of about 6 rows, LIMITATIONS 25; the one-map claim stops %g rows inside '
                                         'the image\'s end, wherever the mask ends it); set with the results in view '
                                         '(regression limits).  Not blind: the start width was fitted on leg 1 (legs 3 '
                                         'and 5 test it in part), leg 7\'s width follows from F (the text\'s 144 um), the '
                                         'gap (17.1 um) and the leg axes were read on this micrograph, so the flanks and '
                                         'the turns are consistency; blind: the scale (the printed centre), the rotation '
                                         '(the image) and V1' % (px_um, le7['margin_rows'])),
                            leg7_claim_end=le7)
    rdm = rd['measured']['side_channels_at_leg_ends']
    dd = [r['image_um'] - r['mask_um'] for r in rdm]
    res['design_side_channels_at_leg_ends'] = dict(
        rows=[dict(r, image_minus_mask_um=r['image_um'] - r['mask_um']) for r in rdm],
        image_minus_mask_um_range=[min(dd), max(dd)],
        note='the design-level reference (sealed start, the caption\'s 140 um) on its one-map registration: reported')
    rec = json.load(open(DEV['aligned']['record']))
    rc = rec['validation']['realised_cutoff']
    res['entrance'] = dict(declared=rc['declared_rows'], worst_over=rc['declared_worst_over'],
                           share_of_F_oversize=rc['declared_share_of_F_oversize'], row1_cutoff_over_D=rc['row1_cutoff_over_D'],
                           note='the start rule\'s rows: computed by the network solve of the drawn device, reported')
    res['physics'] = {k: K.physics_block(json.load(open(DEV[k]['record']))) for k in ('aligned', 'design')}
    res['null_ratio'] = dict(one_map=K.null_ratio(ra1), per_leg={k: K.null_ratio(ra4[k]) for k in LEGS})
    with open(os.path.join(OUT, 'alignment.json'), 'w') as fh:
        json.dump(res, fh, indent=1, default=K.jsonable)
    if not all(c['ok'] for c in checks):
        K.fail('the aligned photomask is not within the alignment limits')
    K.say('measure: the aligned photomask: %d checks within their regression limits (set with the results in view) '
          '-> out/alignment.json' % len(checks))


def _panels(key, T, P6):
    """Return the photomask, its one-map panel, and its per-leg panels, after checking that each claim splits the
    edges as the trace did and that the one-map panel was registered on this mask."""
    d = _need_dxf(key)
    if not (os.path.exists(d['reg_one_map']) and os.path.exists(d['reg_per_leg'])):
        K.fail('run the trace step first')
    r1, r4 = json.load(open(d['reg_one_map'])), json.load(open(d['reg_per_leg']))
    mask = T.Mask(d['dxf'])
    sm1, cl1, pm = split_mask(mask, r1['four_legs'], P6, True, d['preset'])
    n = check_claim(sm1, cl1, pm, r1['four_legs'], 'four_legs')
    le7 = r1['four_legs']['measured']['leg7_claim_end']
    K.say('render (%s): one map: %d / %d edges claimed (cyan) / not (amber: leg 7\'s last %.1f drawn rows, to %d rows '
          'inside the image\'s leg-7 end, and the outlet stretch), as the trace split them' % (
              key, n[0], n[1], le7['rows_before_drawn_end'], le7['margin_rows']))
    sm4, cl4, _ = split_mask(mask, r4['leg5'], P6, False, d['preset'])
    for k in LEGS:
        n = check_claim(sm4, cl4, pm, r4[k], k)
    K.say('render (%s): each leg on its own map: %d / %d edges, as the trace split them' % (key, n[0], n[1]))
    rasters = {K.sha256(IMAGE): K.load_rgb(IMAGE)}
    panels = K.prepare_panels(PANELS, r1, rasters, sm1, claim={'four_legs': cl1})
    for p in panels:
        if p['rec']['mask_sha256'] != mask.sha256:
            K.fail('%s was registered on another mask' % p['key'])
    per_leg = K.prepare_panels(LEGS, r4, rasters, sm4, claim={k: cl4 for k in LEGS})
    return mask, panels, per_leg


def step_render():
    K.say('== render')
    K.use_tree()
    import trace_lib as T
    import panels_gifford6 as P6
    mask, panels, per_leg = _panels('aligned', T, P6)
    le7_margin = panels[0]['rec']['measured']['leg7_claim_end']['margin_rows']
    z1 = K.render_zoom(panels, mask, os.path.join(OUT, 'fig_zoom_overlay'), overlay=True)
    z0 = K.render_zoom(panels, mask, os.path.join(OUT, 'fig_zoom_micrographs'), overlay=False)
    K.say('render: zoomed panel %d x %d px (each image pixel x %d), with and without the overlay' % (
        z1['width_px'], z1['height_px'], z1['magnify']))
    full = K.render_full(mask, panels, os.path.join(OUT, 'fig_full_device'))
    K.say('render: whole device: PDF page %.1f x %.1f in (%g x the device), PNG %d x %d px at %g um/px' % (
        full['pdf']['page_in'][0], full['pdf']['page_in'][1], full['pdf']['magnification'],
        full['png']['size_px'][0], full['png']['size_px'][1], full['png']['um_per_px']))
    stack = K.render_stack(mask, panels, os.path.join(OUT, 'fig_full_device_stacked'), full, rows=STACK_ROWS)
    K.say('render: whole device in %d stacked rows of %s mm (micrograph: %s): PDF page %.1f x %.1f in, '
          'PNG %d x %d px' % (stack['rows'], '/'.join('%.2f' % v for v in stack['row_length_mm']),
                              ', '.join('%s in row %d' % kv for kv in stack['panel_row'].items()),
                              stack['pdf']['page_in'][0], stack['pdf']['page_in'][1], *stack['png']['size_px']))
    mask_d, panels_d, per_leg_d = _panels('design', T, P6)
    zd = K.render_zoom(panels_d, mask_d, os.path.join(DEV['design']['out'], 'fig_zoom_overlay'), overlay=True)
    K.say('render: the reference (design) zoomed panel with its overlay, %d x %d px' % (zd['width_px'], zd['height_px']))
    summary = dict(device=DEV['aligned']['name'], paper='Gifford et al. 2014, Lab Chip 14, 4496-4505, Fig. 6 micrograph',
                   code_stamp=K.STAMP,
                   photomask=dict(file='photomask/%s.dxf' % DEV['aligned']['name'], sha256=mask.sha256),
                   zoom_overlay=z1, zoom_micrographs=z0, full_device=full, full_device_stacked=stack,
                   panels=K.registration_summary(panels), per_leg_registrations=K.registration_summary(per_leg),
                   reference_design=dict(photomask=dict(file='design/photomask/%s.dxf' % DEV['design']['name'],
                                                        sha256=mask_d.sha256), zoom_overlay=zd,
                                         panels=K.registration_summary(panels_d),
                                         per_leg_registrations=K.registration_summary(per_leg_d)),
                   colours=dict(cyan='the mask\'s edges the one-map registration claims (legs 7, 5, 3, 1 and the turns '
                                     'of legs 5, 3, 1)',
                                amber='the edges it does not claim (leg 7\'s drawn rows from %g rows inside the image\'s '
                                      'leg-7 end on, and the outlet stretch)' % le7_margin),
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
    main(sys.argv[1:])
