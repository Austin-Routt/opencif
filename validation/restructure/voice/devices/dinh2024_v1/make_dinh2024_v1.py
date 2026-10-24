"""make_dinh2024_v1.py: device 1 of 4 of the opencif 1.0.0 device kits, device V1 of Dinh et al. 2024 (Lab Chip 24,
913-923), reproduced with opencif. The kit draws two readings of the same device, traces each on Fig. 2B's four
micrographs, and shows each whole.

    python3 make_dinh2024_v1.py            all steps: check, build, trace, render
    python3 make_dinh2024_v1.py render     one step (every run re-checks the code and the inputs first)

The kit builds two photomasks with the same code, the release's opencif under tree/:
  depicted   inputs/dinh2024_v1_depicted.yaml, the printed design (f* 1.04e-4 at the 75 um inlet width), with D 2.17 um
             at every gap except the declared departures. These are gap 1 (the seed), head rows 2-18 (the head as
             Fig. 2A's inset draws it), and the papers' rectilinear transition, which carries reverse flow. The
             design has Fig. 2B's square, rounded quarter turns and its rounded head-gap corners. It is the reference
             reading, and its figures are out/fig_*.
  aligned    inputs/dinh2024_v1_figure_aligned.yaml, the same head, with the linear section at the cutoff that makes
             the drawn device the printed 75.8 mm long (D 3.67 um, fitted; register INV-1, revision). It is a
             two-cutoff reading of the figures, against Dinh's "desired size cut-off throughout the device", and its
             figures are out/figure_aligned/.

The steps run in this order, and every step stops with exit code 2 if a check fails:
  check    the package's opencif against the release (its code stamp against tree/validation/RELEASE.json) and every
           file of the package against MANIFEST.json, the trace code among them; the inputs' sha256 against
           reference/EXPECT.json, which holds RC8's records under the release's code stamp; each generator input
           against its explicit table (INPUTS_DEPICTED, INPUTS_ALIGNED) in both directions; and the aligned preset's
           prescribed head against the depicted preset's planned head, value for value
  build    both photomasks with opencif's generator, each byte-identical to the release's photomask of the same
           preset (which is RC8's rebuild) and with its exit code (that of RC8's record). The step runs the frozen
           physics' own checks of the drawn device (network F = target to the record's Q2 tolerance, no undeclared row
           outside the cutoff band, the declared departures that RC8's record declares, one fluid component, nothing
           lost) and writes its physics block (the plateau, App. B C7's band below it, the reverse flow, the share of
           F within 1% of D, and the model band's junction terms), the rows drawn (head, linear), and the lengths.
  trace    the four Fig. 2B micrographs registered on each photomask by the trace code (the scale from each panel's
           50 um bar, the rotation measured on the image, two independent registrations, the nulls, and a scale
           scan), with every value equal to RC8's trace record (reference/, from provenance/rc7_traces/results/dinh/)
  render   for each reading, after checking that each panel's registration names this photomask by its sha256 and
           that the stacked rows hold every micrograph whole, the zoomed panels with and without the overlay, the
           whole device, and the whole device in stacked rows (PDF + PNG each); then out/registration.json (with each
           reading's physics block), out/run_log.txt (every run's lines), and out/SHA256SUMS

The figures add no text, so that the micrographs keep the publisher's own labels. Every number is in
out/registration.json, and README.md is written from it and from the build records."""
import json
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import devicekit as K   # noqa: E402

RASTERS = os.path.join(HERE, 'inputs')
OUT = os.path.join(HERE, 'out')
REF = os.path.join(HERE, 'reference')
with open(os.path.join(REF, 'EXPECT.json')) as _fh:
    EXPECT = json.load(_fh)
DEV = {
    'depicted': dict(name='dinh2024_v1_depicted', out=OUT, tag='rc7', label='the shipped mask dinh2024_v1_depicted',
                     trace_ref='dinh_results_rc7.json', stack_rows=11),
    'aligned': dict(name='dinh2024_v1_figure_aligned', out=os.path.join(OUT, 'figure_aligned'), tag='fa',
                    label='the figure-aligned mask dinh2024_v1_figure_aligned', trace_ref='dinh_results_fa.json',
                    stack_rows=5),
}
for d in DEV.values():
    e = EXPECT['builds'][d['name']]
    d.update(exit=e['exit_code'], declared=tuple(e['declared']), dxf_sha256=e['dxf_sha256'])
    d['preset'] = os.path.join(RASTERS, d['name'] + '.yaml')
    d['photomask'] = os.path.join(d['out'], 'photomask')
    d['dxf'] = os.path.join(d['photomask'], d['name'] + '.dxf')
    d['record'] = os.path.join(d['photomask'], d['name'] + '.json')
    d['registration'] = os.path.join(d['out'], 'registration_full.json')
    d['ref_registration'] = os.path.join(REF, d['trace_ref'])
FILES = {os.path.join(HERE, k): v for k, v in EXPECT['inputs_sha256'].items()}
# Fig. 2 and Fig. 3 (Fig. 3 is not traced)
IMAGES = [os.path.join(RASTERS, 'dinh5-000.png'), os.path.join(RASTERS, 'dinh6-000.png')]

# The generator inputs, written out. The codes are P, printed in the paper; M, measured on its figures; D, derived from
# printed values through the frozen model (or fitted to a measurement, which the entry says); A, assumed where the
# paper is silent; and G, the generator's own choice. The check step asserts that each preset holds exactly these
# values.
_COMMON = [
    ('targets.filtrate_fraction', 0.5592, 'D', 'outlet identity at the printed outlet widths w_c 52, w_S(out) 44 um (measured F 0.505 +- 0.033, Fig. 4)'),
    ('process.depth_um', 140, 'P', 'Fig. 2 caption: "Depth of all channels was 140 um"'),
    ('process.gap_um', 19, 'P', 'Fig. 2 caption: "All filtration gaps in the linear section were 19 um-wide"'),
    ('process.post_length_um', 26.1, 'M', 'g_p (glossary: post width, along the flow) = row pitch 45.1 +- 0.2 um on Fig. 2B minus G 19'),
    ('process.post_width_um', 54, 'M', 'l_g (link length, across the flow): Fig. 2B post columns 106.06 um apart minus the printed centre 52'),
    ('process.min_feature_um', 19, 'A', 'Dinh print no floor; set to the gap G (Xia 2016\'s rule)'),
    ('centre_channel.taper', 'constant_shear', 'A', 'Dinh print the end widths only'),
    ('centre_channel.width_in_um', 75, 'P', 'Fig. 2 caption: "narrowed from 75 um near the inlet to 52 um near the outlet"'),
    ('centre_channel.width_min_um', 52, 'P', 'ibid.'),
    ('layout.style', 'linear', 'P', 'Fig. 2A, Fig. 3A (~3.2 mm x 75.8 mm)'),
    ('meander.enabled', True, 'P', 'Fig. 2 caption: "the meandering side channel segments"; Fig. 2A'),
    ('meander.style', 'arms', 'M', 'Fig. 2A inset: three goal-post arms, then hairpins'),
    ('meander.arms', 3, 'M', 'Fig. 2A inset'),
    ('meander.channel_um', 28, 'P', 'Fig. 2 caption: "Width of the meandering side channel segments was 28 um"'),
    ('meander.wall_um', 27, 'M', 'Fig. 2A inset'),
    ('meander.head_width_um', 3200, 'P', 'p. 916: footprint ~3.2 mm x 75.8 mm'),
    ('meander.seam_side_width_um', 21, 'P', 'Fig. 2 caption: w_S(in) = 21 um'),
    ('meander.seal', 'allow', 'G', 'draw the depicted head although it does not seal gap 1 at f* (declared: "gap 1 above f*")'),
    ('meander.transition', 'paper', 'M', 'Fig. 2B second panel: the papers\' rectilinear transition (declared departure)'),
    ('meander.rectilinear_rows', 3, 'M', 'Fig. 2B second panel (three wall blocks)'),
    ('meander.seam_at', 'post', 'M', 'Fig. 2B second panel: w_S(in) is the side channel beside the first stadium post'),
    ('meander.turns', 'square', 'M', 'rc7 (DINH-1): Fig. 2B first panel: the head\'s quarter turns are square corners, rounded'),
    ('meander.corner_radius_um', 7.0, 'M', 'the outer corners\' rounding (6.0 / 8.6 um left / right, provenance/rc7_scripts/measure_p1.py)'),
    ('meander.inner_corner_radius_um', 10.0, 'M', 'the inner corners (the walls\' ends): the best corner radius, 9-11 um'),
    ('meander.junction_radius_um', 11.0, 'M', 'rc7 (DINH-7): the head gaps\' entrance corners on the trunk, about 7-15 um '
                                              '(registration-dependent); drawing only (c_ent is an external input)'),
    ('inlet.style', 'dinh_filter', 'M', 'Fig. 2A (two-stage aggregate filter); radii read on Fig. 2A, scale not established'),
    ('inlet.port_radius_um', 1020, 'M', 'Fig. 2A main panel at the arm span taken as 3.2 mm (that scale is not established)'),
    ('inlet.stem_half_width_um', 600, 'M', 'Fig. 2A inset: the filter\'s stem beside which the arms rise'),
    ('outlet.style', 'teardrop', 'M', 'Fig. 2A'),
    ('outlet.port_radius_um', 980, 'M', 'Fig. 2A main panel, as the inlet port (scale not established)'),
    ('outlet.centre_port_radius_um', 956, 'M', 'ibid.'),
    ('outlet.flare_deg', 14.5, 'M', 'Fig. 2B fourth panel, both walls'),
    ('outlet.flare_start_um', 65.4, 'M', 'Fig. 2B fourth panel: straight 65.4 +- 2 um from the last gap\'s centre'),
    ('reproduction.margin', 15.0, 'G', 'the head\'s lateral margin inside the 3.2 mm width'),
    ('reproduction.R0_arms', 1281.0, 'M', 'Fig. 2A inset: arms at the chip edge'),
    ('reproduction.ws_rect', 28.0, 'M', 'Fig. 2B second panel: the 28 um side channel beside the wall blocks'),
    ('gates.mode', 'reproduce', 'G', 'draw what the paper depicts and report the envelope (strict mode would refuse: below the kappa plateau)'),
    ('operation.flow_uL_min', 30, 'P', 'p. 915: 30 uL/min'),
]
INPUTS_DEPICTED = _COMMON + [
    ('targets.cutoff_um', 2.1741549811760246, 'D', 'printed f*gap 1.04e-4 at the inlet width 75 um, through kappa_g (App. B)'),
    ('reproduction.depicted_head', 'dinh', 'M', 'Fig. 2A inset unit lengths (18 rows, 3 arms), derived for the square-turn head'),
]
INPUTS_ALIGNED = _COMMON + [
    ('targets.cutoff_um', 3.67, 'D', 'fitted: the drawn device (head, linear section, inlet with its filter, outlets) is '
                                     'the printed 75.8 mm long (75.71 mm drawn, 874 linear rows; register INV-1, revision); '
                                     'Dinh\'s 4.7 um beads imply 2.3-3.9 um on Gifford\'s soft edge; the figures\' row '
                                     'counts give 4.28-4.79 um, their proportion given up'),
    ('reproduction.NM_force', 219, 'D', 'the depicted preset\'s planned head rows'),
    ('reproduction.f_head', None, 'D', 'the depicted preset\'s planned head, value for value (checked against its plan)'),
]
PANELS = ['p1', 'p2', 'p3', 'p4']      # Fig. 2B, top to bottom: head, transition, linear section, outlet


def depicted_head():
    """Return the depicted preset's planned head fractions and its number of head rows, which the aligned preset must
    prescribe exactly."""
    K.use_tree()
    from opencif import spec as SP
    from opencif.physics import schedule as S
    sp, _ = SP.to_engine(SP.load(DEV['depicted']['preset']))
    S.clear_caches()
    sc = S.plan(sp)
    return [float(v) for v in sc.f[:int(sc.NM)]], int(sc.NM)


def step_check():
    K.say('== check')
    K.check_code()
    K.say('reference/EXPECT.json (the frozen tree\'s records): sha256 %s' % K.sha256(os.path.join(REF, 'EXPECT.json'))[:16])
    if EXPECT['code_stamp'] != K.STAMP:
        K.fail('reference/EXPECT.json is for code %s, the package for %s' % (EXPECT['code_stamp'], K.STAMP))
    K.check_inputs(FILES)
    K.say('generator input (depicted, the printed reading):')
    K.check_preset(DEV['depicted']['preset'], INPUTS_DEPICTED)
    f_head, nm = depicted_head()
    table = [(p_, (f_head if p_ == 'reproduction.f_head' else v_), c_, s_) for p_, v_, c_, s_ in INPUTS_ALIGNED]
    nm_table = dict((p_, v_) for p_, v_, _, _ in INPUTS_ALIGNED)['reproduction.NM_force']
    nm_ref = EXPECT['builds'][DEV['depicted']['name']]['NM']
    if not nm == nm_table == nm_ref:
        K.fail('the depicted preset plans %d head rows; the aligned table says %d, the tree\'s record %d' % (nm, nm_table, nm_ref))
    K.say('generator input (figure-aligned; its head is the depicted plan\'s %d rows):' % nm)
    K.check_preset(DEV['aligned']['preset'], table)


def step_build():
    K.say('== build')
    for key in ('depicted', 'aligned'):
        d = DEV[key]
        K.say('build (%s):' % key)
        _, _, rec = K.build(d['preset'], d['photomask'], d['dxf_sha256'], d['exit'])
        dp = rec['validation'].get('drawn_parameters', {})
        K.say('build (%s): N %d gap rows: %d in the head, %d linear; F %.5f; footprint %s mm; array %s mm, head %s mm' % (
            key, rec['design']['N'], rec['design']['NM'], rec['design']['N'] - rec['design']['NM'], rec['design']['F'],
            dp.get('footprint_mm'), dp.get('array_length_mm'), dp.get('head_length_mm')))
        # F is checked to the record's own Q2 tolerance (1e-3), because the head is drawn as depicted (prescribed row
        # by row) rather than closed on F.
        K.check_build_physics(rec, d['declared'], f_tol=1e-3)


def _need_dxf(key):
    """Return the photomask's entry, and stop unless out/ holds its photomask with the reference sha256."""
    d = DEV[key]
    if not os.path.exists(d['dxf']) or K.sha256(d['dxf']) != d['dxf_sha256']:
        K.fail('out/ holds no %s photomask matching its reference: run the build step' % key)
    return d


def step_trace():
    K.say('== trace')
    K.use_tree()
    import panels_dinh as PD
    import panels_dinh_v2 as PV
    for key in ('depicted', 'aligned'):
        d = _need_dxf(key)
        work = os.path.join(OUT, 'trace_work')            # the trace code's own figures (captioned raster crops),
        os.makedirs(work, exist_ok=True)                  # which are removed afterwards, whatever happens
        try:
            ctx = PD.setup(K.TREE, RASTERS, preset=d['preset'], dxf=d['dxf'], record=d['record'], label=d['label'], tag=d['tag'])
            res, t0 = {}, time.time()
            res['p4'] = PV.panel4(ctx, work, True)
            res['p3'] = PV.panel3(ctx, work, res['p4'], True)
            res['p2'] = PV.panel2(ctx, work, True)
            res['p1'] = PV.panel1(ctx, work, True)
            K.say('trace (%s): Fig. 2B panels 4, 3, 2, 1 registered, %.0f s' % (key, time.time() - t0))
            K.check_trace(res, d['ref_registration'], d['registration'])
        finally:
            shutil.rmtree(work, ignore_errors=True)


def _render(key, T):
    """Draw one reading's figures, after checking that each panel was registered on its photomask, and return the
    reading's entry of out/registration.json."""
    d = _need_dxf(key)
    if not os.path.exists(d['registration']):
        K.fail('run the trace step first')
    res = json.load(open(d['registration']))
    mask = T.Mask(d['dxf'])
    rasters = {K.sha256(p): K.load_rgb(p) for p in IMAGES}
    panels = K.prepare_panels(PANELS, res, rasters, mask)
    for p in panels:
        if p['rec']['mask_sha256'] != mask.sha256:
            K.fail('%s was registered on another mask' % p['key'])
    z1 = K.render_zoom(panels, mask, os.path.join(d['out'], 'fig_zoom_overlay'), overlay=True)
    z0 = K.render_zoom(panels, mask, os.path.join(d['out'], 'fig_zoom_micrographs'), overlay=False)
    K.say('render (%s): zoomed panels %d x %d px (each image pixel x %d), with and without the overlay' % (
        key, z1['width_px'], z1['height_px'], z1['magnify']))
    full = K.render_full(mask, panels, os.path.join(d['out'], 'fig_full_device'))
    K.say('render (%s): whole device: PDF page %.1f x %.1f in (%g x the device), PNG %d x %d px at %g um/px' % (
        key, full['pdf']['page_in'][0], full['pdf']['page_in'][1], full['pdf']['magnification'],
        full['png']['size_px'][0], full['png']['size_px'][1], full['png']['um_per_px']))
    stack = K.render_stack(mask, panels, os.path.join(d['out'], 'fig_full_device_stacked'), full, rows=d['stack_rows'])
    K.say('render (%s): whole device in %d stacked rows of %s mm (micrographs: %s): PDF page %.1f x %.1f in, '
          'PNG %d x %d px' % (key, stack['rows'], '/'.join('%.2f' % v for v in stack['row_length_mm']),
                              ', '.join('%s in row %d' % kv for kv in stack['panel_row'].items()),
                              stack['pdf']['page_in'][0], stack['pdf']['page_in'][1], *stack['png']['size_px']))
    summ = K.registration_summary(panels)
    for k, v in summ.items():
        if v['below_criterion']:
            K.say('render (%s): panel %s: null ratio below %.1f for %s (%s): flagged' % (
                key, k, K.NULL_RATIO_CRITERION, ', '.join(v['below_criterion']),
                ', '.join('%s %.2f' % kv for kv in sorted(v['null_ratio'].items()))))
    return dict(device=d['name'], photomask=dict(file=os.path.relpath(d['dxf'], OUT), sha256=mask.sha256),
                zoom_overlay=z1, zoom_micrographs=z0, full_device=full, full_device_stacked=stack,
                panels=summ, physics=K.physics_block(json.load(open(d['record']))))


def step_render():
    K.say('== render')
    K.use_tree()
    import trace_lib as T
    summary = dict(paper='Dinh et al. 2024, Lab Chip 24, 913-923, Fig. 2B', code_stamp=K.STAMP,
                   depicted=_render('depicted', T), figure_aligned=_render('aligned', T),
                   coordinates='mask um: x along the device, y across; image px: (column, row) of the raster; '
                               'image px = origin + x ex + y ey (map_chamfer)')
    with open(os.path.join(OUT, 'registration.json'), 'w') as fh:
        json.dump(summary, fh, indent=1)


STEPS = dict(check=step_check, build=step_build, trace=step_trace, render=step_render)


def main(argv):
    todo = argv or list(STEPS)
    os.makedirs(OUT, exist_ok=True)
    try:
        if 'check' not in todo:
            step_check()                               # every run re-checks the code and the inputs
        for s in todo:
            STEPS[s]()
        K.say('done: %s' % ', '.join(todo))
    except Exception as e:                             # a check's failure is logged by K.fail; anything else here
        K.say('STOPPED: %s: %s' % (type(e).__name__, e))
        raise
    finally:
        with open(os.path.join(OUT, 'run_log.txt'), 'a') as f:      # every run's lines kept (package check, 28)
            f.write('\n'.join(K.LOG) + '\n')
        if 'render' in todo:
            K.write_sums(OUT)


if __name__ == '__main__':
    main(sys.argv[1:])
