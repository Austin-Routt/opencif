"""make_gifford2014_fig4.py: device 2 of 4 of the opencif 1.0.0 device kits, the Fig. 4 bead device of Gifford et al.
2014 (Lab Chip 14, 4496-4505), reproduced with opencif and aligned with the Fig. 4 inset by the nested reading.

The nested reading is D1 of the adjudication of the rc7 reviews (2026-10-20; register GIF4-1 / INV-2). The legs are
nested, as the caption ("the second and fourth legs") and the schematic say, as the inset's two co-directional arrows
suggest, and as the inset's own taper direction confirms (V13). The arrows are the weakest line of this evidence,
because they are placed symmetrically about leg 4 and so may both mark leg 4. The device starts by Gifford's start
rule at an 8 um first side channel (rows 1-3 above D, row 1 at 2.0 D, 0.5% of F, declared), and its gap is the one
measured on the inset. One feature is not reproduced, and the kit states so. Leg 2's side channels beside the inset
read about 18 um per side wider on the image than the printed single-f schedule draws them (register GIF4-1, open).

    python3 make_gifford2014_fig4.py            all steps: check, build, trace, taper, measure, render
    python3 make_gifford2014_fig4.py render     one step (every run re-checks the code and the inputs first)

The kit builds five photomasks with the same code, the release's opencif under tree/:
  aligned    inputs/gifford2014_fig4_aligned.yaml, Gifford's printed design with every printed value, his start rule
             at the 10 : 1 aspect ratio that he states (an 8 um first side channel), the four legs nested, leg 4 ending
             2 mm short of the others as the inset's leads place it, and the gap measured on the inset. It is the
             package's device, and its figures are out/fig_*.
  design     inputs/gifford2014_fig4_design.yaml, the design-level preset (sealed start, nested legs, 20.5 um gap),
             which is the reference (out/design/).
  pathorder  inputs/gifford2014_fig4_pathorder.yaml, the rejected reading B, which was rc7's aligned preset. Its legs
             lie side by side in path order, so that the inset's neighbor would be leg 3. Its widths fit. However, the
             neighbor's taper runs against leg 3's flow (V13 FAIL), and the caption, the schematic, and the arrows say
             nested. It is built, traced, and shown as a figure only (out/rejected_pathorder/).
  start51    inputs/gifford2014_fig4_start51.yaml, the rejected reading C, with nested legs and a 51 um start fitted to
             leg 2's widths. The inset's beads rule it out (register INV-2), so it is built, traced, and shown as a
             figure only (out/rejected_start51/).
  fold20     inputs/gifford2014_fig4_fold20.yaml, the 20 um folded-head variant (A. H. Routt, 2026-10-19, who asked
             to keep the 8 um start and to show the folded head). It draws spec v5's FOLD start at a 20 um floor (the
             lab's later guideline, Lezzar et al. 2022), with every gap at D, legs 1-2 at the generator's own spacing
             for the head, and legs 2-4 on reading B's rows, so it is a path-order figure. It is not Gifford's device,
             because he drew no head and the inset's gaps read 17.2 um. It is built, traced, and shown as a figure only
             (out/variant_fold20/).
The three variants (pathorder, start51, and fold20) take their references from RC8's builds of them.

The steps run in this order, and every step stops with exit code 2 if a check fails:
  check    the package's opencif against the release (its code stamp against tree/validation/RELEASE.json) and every
           file of the package against MANIFEST.json, the trace code among them; the inputs' sha256 against
           reference/EXPECT.json, which holds RC8's records under the release's code stamp; and each generator input
           against its explicit table (INPUTS_*) in both directions
  build    the five photomasks, built with opencif's generator, each byte-identical to the release's photomask of the
           same spec (which is RC8's) and with its exit code. For each, the step runs the frozen physics' own checks
           of the drawn device (network F = target, no undeclared row outside the cutoff band, one fluid component,
           nothing lost) and writes its physics block (the share of F within 1% of D, the declared departures, and
           the model band's junction terms). The aligned and the rejected devices declare their start rows
           ('floor start').
  trace    the inset registered on each photomask by the trace code, with every value equal to RC8's trace record
           (reference/, from provenance/rc7_traces/results/gifford/; only the aligned and design records are anchors)
  taper    V13, the taper-direction sign test (tree/provenance/rc7_traces/code/taper_sign.py), on each photomask's
           neighboring leg (wall window +-30 um) and on leg 4's array part (+-12 um, the control), from each trace's
           own map. Every result must equal RC8's (reference/v13_*.json), the aligned reading's neighbor must PASS,
           reading B's must FAIL, and every leg-4 control must PASS -> out/taper_direction.json
  measure  on each photomask, the widths that the inset shows against the drawn ones (leg 4 and the leg beside it),
           the neighbor's walls (from the trace record; a wall with fewer than 20 samples is listed but not checked),
           and leg 4's posts' flanks along the flow at the trace's own search window. The limits hold on the aligned
           device's claim (leg 4) only. They are one image pixel (the record's scale) for leg 4's side channels,
           0.5 um for the flanks' median, and 1 um for the claim's median, and they are regression limits, set with
           the results in view. The neighbor's widths and walls are reported for every photomask, and for the aligned
           device they are stated as not reproduced (GIF4-1) -> out/alignment.json
  render   the aligned device's figures (out/fig_zoom_overlay, fig_zoom_micrographs, fig_full_device,
           fig_full_device_stacked; PDF + PNG each), the reference's and the three other photomasks' zoom overlays
           (out/design/, out/rejected_pathorder/, out/rejected_start51/, out/variant_fold20/), the V13 figure
           (out/fig_taper_direction, widths only and no paper raster), out/registration.json, out/run_log.txt (every
           run's lines appended), and out/SHA256SUMS

The figures add no text, so that the micrograph keeps the publisher's own labels. Every number is in
out/registration.json, out/alignment.json, and out/taper_direction.json, and README.md is written from them and from
the build records."""
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
IMAGE = os.path.join(RASTERS, 'gif5-001.png')
with open(os.path.join(REF, 'EXPECT.json')) as _fh:
    EXPECT = json.load(_fh)
DEV = {
    'aligned': dict(name='gifford2014_fig4_aligned', out=OUT, tag='aligned',
                    label='the aligned mask gifford2014_fig4_aligned (nested legs, RC8)',
                    trace_ref='gifford_results_aligned.json'),
    'design': dict(name='gifford2014_fig4_design', out=os.path.join(OUT, 'design'), tag='rc7',
                   label='the shipped mask gifford2014_fig4_design', trace_ref='gifford_results_rc7.json'),
    'pathorder': dict(name='gifford2014_fig4_pathorder', out=os.path.join(OUT, 'rejected_pathorder'), tag='pathorder',
                      label="the rejected variant gifford2014_fig4_pathorder (rc7's reading B, legs in path order)",
                      trace_ref='gifford_results_pathorder.json'),
    'start51': dict(name='gifford2014_fig4_start51', out=os.path.join(OUT, 'rejected_start51'), tag='start51',
                    label='the rejected variant gifford2014_fig4_start51 (51 um start, nested)',
                    trace_ref='gifford_results_start51.json'),
    'fold20': dict(name='gifford2014_fig4_fold20', out=os.path.join(OUT, 'variant_fold20'), tag='fold20',
                   label='the folded-head variant gifford2014_fig4_fold20 (20 um floor, legs in path order)',
                   trace_ref='gifford_results_fold20.json'),
}
ORDER = ('design', 'aligned', 'pathorder', 'start51', 'fold20')
for d in DEV.values():
    e = EXPECT['builds'][d['name']]
    d.update(exit=e['exit_code'], declared=tuple(e['declared']), dxf_sha256=e['dxf_sha256'])
    d['preset'] = os.path.join(RASTERS, d['name'] + '.yaml')
    with open(d['preset']) as _fh:
        _p = K.yaml_load(_fh)['process']
    d['post_max'] = (float(_p['post_length_um']) + 1.0, float(_p['post_width_um']) + 1.0)   # the trace's post filter
    d['photomask'] = os.path.join(d['out'], 'photomask')
    d['dxf'] = os.path.join(d['photomask'], d['name'] + '.dxf')
    d['record'] = os.path.join(d['photomask'], d['name'] + '.json')
    d['registration'] = os.path.join(d['out'], 'registration_full.json')
    d['ref_registration'] = os.path.join(REF, d['trace_ref'])
    d['v13'] = {w: os.path.join(REF, 'v13_%s_%s.json' % (d['tag'], w)) for w in ('neighbour', 'leg4')}
FILES = {os.path.join(HERE, k): v for k, v in EXPECT['inputs_sha256'].items()}

# The generator inputs, written out. The codes are P, printed in the paper; M, measured on its figures; D, derived from
# printed values through the frozen model (or fitted to a measurement, which the entry says); A, assumed where the
# paper is silent; and G, the generator's own choice.
_COMMON = [
    ('targets.cutoff_um', 5.203174831309222, 'D', 'printed f_gap 3.34e-4 at w_c 100 um (Fig. 4 caption), through kappa_g (App. B)'),
    ('targets.filtrate_fraction', 0.90644, 'D', 'outlet identity at the printed final side width 300 um (Fig. 4 caption)'),
    ('process.depth_um', 80, 'P', 'Fig. 4 caption: "The depth of the device was 80 um +- 3 um"'),
    ('process.post_width_um', 59.5, 'M', 'l_g (across the flow): each leg\'s post columns 159.5 um apart minus the printed 100 um centre'),
    ('centre_channel.taper', 'constant', 'P', 'Gifford\'s design: a constant central channel (p. 4498, Fig. 2)'),
    ('centre_channel.width_in_um', 100, 'P', 'Fig. 4 caption: "the width of the central channel was 100 um"'),
    ('layout.style', 'paperclip', 'M', 'Fig. 4 schematic (not drawn to scale, p. 4502): four legs, the ports at the top'),
    ('layout.legs', 4, 'P', 'Fig. 4 caption ("the second and fourth legs of the device"); the schematic\'s four legs'),
    ('layout.leg_spacing_um', 818, 'M', 'Fig. 4 inset: the centre lines of leg 4 and the leg beside it 817.9 um apart'),
    ('meander.enabled', False, 'P', 'Gifford\'s devices have no meander'),
    ('inlet.style', 'plain', 'M', 'Fig. 4 schematic: a plain round port'),
    ('inlet.port_radius_um', 300, 'A', 'the schematic is not to scale'),
    ('outlet.style', 'separate', 'M', 'Fig. 4 schematic: the centre and side outputs at separate ports'),
    ('outlet.port_radius_um', 750, 'A', 'the schematic is not to scale'),
    ('outlet.centre_port_radius_um', 300, 'A', 'ibid.'),
    ('gates.mode', 'reproduce', 'G', 'draw what the paper depicts and report what the process limits flag'),
    ('reproduction.stem_half', 300.0, 'G', 'an inlet-filter parameter (the replica record\'s); unused without a filter'),
    ('reproduction.W', 2400.0, 'G', 'the head envelope\'s width (the replica record\'s); this device has no head'),
]
INPUTS_DESIGN = _COMMON + [
    ('process.gap_um', 20.5, 'M', 'the record\'s measurement on Fig. 4 (p. 4496: the gap size "~20 um"); GIF4-3'),
    ('process.post_length_um', 25.3, 'M', 'g_p (along the flow) = the inset\'s row pitch 45.82 um minus G 20.5'),
    ('process.min_feature_um', 4.5, 'G', 'the design preset\'s floor (chosen so its sealed start\'s 4.86 um first channel passes)'),
    ('process.min_radius_um', 1.8, 'G', 'ibid.'),
    ('process.fillet_um', 2.0, 'G', 'fillet over the inlet and outlet zones'),
    ('reproduction.pc_last_short', 2000.0, 'M', 'Fig. 4 inset: leg 4\'s post array ends before the neighbour\'s (reproduction-only)'),
    ('reproduction.wall_M', 4.5, 'G', 'a meander-head parameter (the replica record\'s); this device has no head'),
]
_GIFFORD_START = [
    ('process.gap_um', 17.2, 'M', 'the inset\'s posts\' flanks along the flow at half-windows of 3.5 and 6 px (GIF4-2, GIF4-3); '
                                  'checked by this package\'s measure step at the trace window (at most half the gap)'),
    ('process.post_length_um', 28.6, 'M', 'g_p (along the flow) = the inset\'s row pitch 45.8 um minus G 17.2'),
    ('process.min_feature_um', 8, 'A', 'the floor at the start width (Gifford\'s "10 : 1 aspect ratio limitation", p. 4503, at 80 um)'),
    ('process.fillet_um', 2.5, 'G', 'fillet over the inlet and outlet zones (below half the 8 um floor)'),
    ('start.style', 'floor', 'P', 'Gifford\'s start rule (p. 4498): w_s(1) "close to or equal to zero, depending on the device '
                                  'production capabilities available"'),
    ('start.side_width_um', 8, 'A', 'his "10 : 1 aspect ratio limitation" (p. 4503) at the 80 um depth: an inference'),
]
INPUTS_ALIGNED = _COMMON + _GIFFORD_START + [
    ('reproduction.pc_last_short', 2000.0, 'M', 'Fig. 4 inset: leg 4\'s post array ends before the neighbour\'s (nested layout, '
                                                'as the design preset; reproduction-only)'),
]
INPUTS_PATHORDER = _COMMON + _GIFFORD_START + [
    ('layout.order', 'path', 'A', 'the REJECTED reading B (register GIF4-1, INV-2): the legs side by side in path order, so leg 3 '
                                  'would be beside leg 4 on the inset; rejected by the inset\'s taper direction (V13), the '
                                  'caption\'s "second and fourth", the schematic\'s nesting (and the arrows, read as one per leg)'),
    ('reproduction.pc_last_short', 8200.0, 'D', 'fitted: leg 4\'s array ends 8.2 mm before the others, where the row-fraction '
                                                'law puts the neighbour\'s measured 124.7 um side channels (register INV-2)'),
]
INPUTS_START51 = _COMMON + [
    ('process.gap_um', 17.2, 'M', 'as the aligned preset (the inset\'s posts\' flanks)'),
    ('process.post_length_um', 28.6, 'M', 'as the aligned preset'),
    ('process.min_feature_um', 4.5, 'G', 'the design preset\'s floor (unchanged; the side channels start at 51 um)'),
    ('process.min_radius_um', 1.8, 'G', 'ibid.'),
    ('process.fillet_um', 2.0, 'G', 'fillet over the inlet and outlet zones'),
    ('start.style', 'floor', 'P', 'Gifford\'s start rule (p. 4498)'),
    ('start.side_width_um', 51.0, 'D', 'fitted: the start that puts leg 2\'s side channels at the inset\'s 124.7 um (M) through the '
                                      'row-fraction law, legs nested; REJECTED by the inset\'s beads (register INV-2)'),
    ('reproduction.pc_last_short', 2000.0, 'M', 'Fig. 4 inset: leg 4\'s post array ends before the neighbour\'s (as the design preset)'),
    ('reproduction.wall_M', 4.5, 'G', 'a meander-head parameter (the replica record\'s); this device has no head'),
]
INPUTS_FOLD20 = [
    ('targets.cutoff_um', 5.203174831309222, 'D', 'printed f_gap 3.34e-4 at w_c 100 um (Fig. 4 caption), through kappa_g (App. B)'),
    ('targets.filtrate_fraction', 0.90644, 'D', 'outlet identity at the printed final side width 300 um (Fig. 4 caption)'),
    ('process.depth_um', 80, 'P', 'Fig. 4 caption: "The depth of the device was 80 um +- 3 um"'),
    ('process.gap_um', 20.5, 'M', 'the design preset\'s (the record\'s measurement on Fig. 4); the inset\'s 17.2 um is below this '
                                  'variant\'s 20 um floor'),
    ('process.post_length_um', 25.3, 'M', 'g_p (along the flow) = the inset\'s row pitch 45.82 um minus G 20.5'),
    ('process.post_width_um', 59.5, 'M', 'l_g (across the flow): each leg\'s post columns 159.5 um apart minus the printed 100 um centre'),
    ('process.min_feature_um', 20, 'A', 'the Shevkoplyas lab\'s later guideline (Lezzar et al. 2022: features "less than 20 um are '
                                       'very difficult to mass produce"; max(20, 80/7) at 80 um): not Gifford\'s floor'),
    ('process.fillet_um', 5.0, 'G', 'the corner rule at a 20 um floor (as gifford2014_fig4_as_fabricated_20um)'),
    ('centre_channel.taper', 'constant', 'P', 'Gifford\'s design: a constant central channel (p. 4498, Fig. 2)'),
    ('centre_channel.width_in_um', 100, 'P', 'Fig. 4 caption: "the width of the central channel was 100 um"'),
    ('layout.style', 'paperclip', 'M', 'Fig. 4 schematic (not drawn to scale, p. 4502): four legs, the ports at the top'),
    ('layout.legs', 4, 'P', 'Fig. 4 caption ("the second and fourth legs of the device"); the schematic\'s four legs'),
    ('layout.order', 'path', 'A', 'as the rejected reading B (gifford2014_fig4_pathorder; register INV-2): a path-order figure'),
    ('meander.enabled', True, 'G', 'the variant: spec v5\'s FOLD start, the rows below the floor drawn as head units (Gifford '
                                   'draws no head)'),
    ('meander.style', 'pack', 'G', 'the generator\'s packing head'),
    ('meander.seal', 'require', 'G', 'gap 1 at f* (sigma0 = 0): every gap at D'),
    ('meander.seam_side_width_um', 20, 'A', 'the first straight side channel at the 20 um floor'),
    ('inlet.style', 'plain', 'M', 'Fig. 4 schematic: a plain round port'),
    ('inlet.port_radius_um', 300, 'A', 'the schematic is not to scale'),
    ('outlet.style', 'separate', 'M', 'Fig. 4 schematic: the centre and side outputs at separate ports'),
    ('outlet.port_radius_um', 750, 'A', 'the schematic is not to scale'),
    ('outlet.centre_port_radius_um', 300, 'A', 'ibid.'),
    ('gates.mode', 'reproduce', 'G', 'reproduction-only layout values (leg axes, inlet offset, leg 4\'s early end)'),
    ('reproduction.pc_last_short', 8200.0, 'D', 'as reading B (fitted there to the neighbour\'s 124.7 um side channels)'),
    ('reproduction.stem_half', 300.0, 'G', 'an inlet-filter parameter (the replica record\'s); unused without a filter'),
    ('reproduction.W', 2400.0, 'G', 'the head envelope\'s width (the replica record\'s): the packing head folds within it'),
    ('reproduction.leg_Y', [0.0, 1672.0, 2490.0, 3308.0], 'G', 'legs 1-2 at the generator\'s own paperclip rule for this head '
                                                             '(max(2 hw_max, hw_head + hw_max) + 150 um); legs 2-3 and 3-4 at the '
                                                             'inset\'s 818 um (M); at 818 um the head lies across leg 2 and rc7 '
                                                             'refuses it (register BUILD-2)'),
    ('reproduction.pc_inlet_offset', -12072.7, 'G', 'the legs\' common end 2836.6 um downstream of gap 1: legs 2-4 on reading '
                                                    'B\'s rows'),
]
INPUTS = dict(aligned=INPUTS_ALIGNED, design=INPUTS_DESIGN, pathorder=INPUTS_PATHORDER, start51=INPUTS_START51,
              fold20=INPUTS_FOLD20)
PANELS = ['inset']        # the Fig. 4 inset, the end of leg 4 (the claim) beside the neighboring leg
STACK_ROWS = 3


def leg4_claim(rec):
    """Return the trace's claim, leg 4's edges, as the edges whose mean y lies above the legs' mid line. The trace
    checks this rule against its own claim on every edge (measured.claim_split.disagree_with_mid_line = 0)."""
    mid = rec['measured']['claim_split']['mid_line_um']
    return lambda e: float(K.np.mean(K.np.asarray(e)[:, 1])) > mid


def step_check():
    K.say('== check')
    K.check_code()
    K.say('reference/EXPECT.json (the frozen tree\'s records): sha256 %s' % K.sha256(os.path.join(REF, 'EXPECT.json'))[:16])
    if EXPECT['code_stamp'] != K.STAMP:
        K.fail('reference/EXPECT.json is for code %s, the package for %s' % (EXPECT['code_stamp'], K.STAMP))
    K.check_inputs(FILES)
    for key, what in (('aligned', 'the package\'s device'), ('design', 'the reference'), ('pathorder', 'the rejected reading B'),
                      ('start51', 'the rejected reading C'), ('fold20', 'the folded-head variant')):
        K.say('generator input (%s, %s):' % (key, what))
        K.check_preset(DEV[key]['preset'], INPUTS[key])


def step_build():
    K.say('== build')
    for key in ORDER:
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
    import panels_gifford as PG
    for key in ORDER:
        d = _need_dxf(key)
        work = os.path.join(OUT, 'trace_work')            # the trace code's own captioned figures (raster crops),
        os.makedirs(work, exist_ok=True)                  # which are removed afterwards, whatever happens
        try:
            t0 = time.time()
            rec = PG.run(K.TREE, RASTERS, work, dxf=d['dxf'], record=d['record'], preset=d['preset'],
                         tag=d['tag'], label=d['label'], scan=True)
            K.say('trace (%s): Fig. 4 inset (leg 4, beside leg %s) registered, %.0f s' % (
                key, rec['measured'].get('neighbour_leg'), time.time() - t0))
            K.check_trace(rec, d['ref_registration'], d['registration'])
        finally:
            shutil.rmtree(work, ignore_errors=True)


TAPER_RUNS = (   # (device, leg, window um, x range) of V13, as provenance/rc7_traces/code/run_all_rc7.sh runs it
    [(k, 'neighbour', 30, None) for k in ('aligned', 'design', 'pathorder', 'start51', 'fold20')]
    + [(k, 'leg4', 12, '@array') for k in ('aligned', 'design', 'pathorder', 'start51', 'fold20')])
TAPER_REQUIRED = {('aligned', 'neighbour'): 'PASS', ('pathorder', 'neighbour'): 'FAIL',
                  ('aligned', 'leg4'): 'PASS', ('design', 'leg4'): 'PASS', ('pathorder', 'leg4'): 'PASS',
                  ('start51', 'leg4'): 'PASS', ('fold20', 'leg4'): 'PASS'}


def step_taper():
    """Run V13, which compares the sign of the image's outer wall-to-wall width slope along a leg with the mask's, from
    each trace's own map (the registration record that this package has just made and checked). Every value is
    compared with RC8's record (reference/v13_*.json)."""
    K.say('== taper')
    import subprocess
    script = os.path.join(K.TRACE, 'taper_sign.py')
    work = os.path.join(OUT, 'taper_work')
    os.makedirs(work, exist_ok=True)
    out = dict(test='V13, the taper-direction sign test (taper_sign.py): PASS when the image\'s width slope has the mask\'s '
                    'sign beyond 3 SE with coverage >= 0.5, >= 30 positions and a mask slope of at least 1e-4 um/um; FAIL '
                    'with the opposite sign beyond 3 SE; else INCONCLUSIVE.  The SE is the larger of the Newey-West SE '
                    '(Bartlett kernel, lag 100 um, beyond two row pitches) and the independent-sample SE: neighbouring '
                    'positions share pixels, so their residuals are correlated and the independent-sample SE alone is too '
                    'small.  The magnitude ratio is reported, never judged '
                    '(the image\'s slope moves by ~15 % with the sample grid).  The wall-search window must exceed the '
                    'largest wall offset: +-30 um on the neighbour (the aligned reading\'s walls lie {RANGE} um outside '
                    'the mask\'s, as this run measured them), +-12 um on leg 4\'s array part.',
               runs={})
    ok = True
    try:
        for key, leg, win, xr in TAPER_RUNS:
            d = _need_dxf(key)
            if not os.path.exists(d['registration']):
                K.fail('run the trace step first')
            res = os.path.join(work, 'v13_%s_%s.json' % (d['tag'], leg))
            cmd = [sys.executable, '-I', script, IMAGE, d['registration'], d['dxf'], leg, '%s_%s' % (d['tag'], leg), res,
                   '%g' % win, '5', '3'] + ([xr] if xr else [])
            p = subprocess.run(cmd, capture_output=True, text=True)
            if p.returncode or not os.path.exists(res):
                K.fail('V13 %s %s did not run: %s' % (key, leg, p.stderr[-300:]))
            r = json.load(open(res))
            ref = json.load(open(d['v13'][leg]))
            same = (r['verdict'] == ref['verdict'] and r['positions'] == ref['positions'] and r['eligible'] == ref['eligible']
                    and all(abs((r.get(k) or 0) - (ref.get(k) or 0)) <= 1e-9 for k in ('slope_image', 'slope_mask', 'se_image',
                                                                                         'se_image_iid', 'se_image_hac')))
            req = TAPER_REQUIRED.get((key, leg))
            good = same and (req is None or r['verdict'] == req)
            ok = ok and good
            K.say('taper (%s, %s, +-%g um): %s, %d of %d eligible positions, image %s, mask %s%s; %s the frozen tree\'s' % (
                key, leg, win, r['verdict'], r['positions'], r['eligible'],
                ('%+.5f +- %.5f (%.1f SE)' % (r['slope_image'], r['se_image'], r['z_image'])) if 'slope_image' in r else '-',
                ('%+.5f' % r['slope_mask']) if 'slope_mask' in r else '-', (' (required %s)' % req) if req else '',
                'equal to' if same else 'NOT EQUAL TO'))
            if (key, leg) == ('aligned', 'neighbour') and r.get('rows'):
                a_ = [abs(x[k_]) for x in r['rows'] for k_ in ('off_low', 'off_high')]
                out['test'] = out['test'].replace('{RANGE}', '%.1f-%.1f' % (min(a_), max(a_)))
            out['runs']['%s_%s' % (key, leg)] = dict((k, v) for k, v in r.items() if k != 'rows')
            out['runs']['%s_%s' % (key, leg)].update(required=req, equal_to_reference=same)
            shutil.copy(res, os.path.join(DEV[key]['out'], 'taper_%s.json' % leg))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    out['test'] = out['test'].replace('{RANGE}', 'about 18')
    with open(os.path.join(OUT, 'taper_direction.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    if not ok:
        K.fail('V13: a result differs from the frozen tree\'s, or a required verdict is not met')
    K.say('taper: every V13 result equal to the frozen tree\'s; the aligned reading\'s neighbour PASS, reading B\'s FAIL, '
          'every leg-4 control PASS -> out/taper_direction.json')


def step_measure():
    K.say('== measure')
    K.use_tree()
    import trace_lib as T
    import gifford_alignment as GA
    raster = T.Raster(IMAGE, EXPECT['inputs_sha256']['inputs/gif5-001.png'][:16])
    res = dict(reading=('the nested reading (D1 of the adjudication of the rc7 reviews; register GIF4-1 / INV-2): the printed '
                        'design with every printed value, the legs nested as the caption ("the second and fourth legs"), the '
                        'schematic and Fig. 6 say, the inset\'s co-directional arrows suggest (the weakest line) and the inset\'s taper direction confirms '
                        '(V13), Gifford\'s start rule at an 8 um first side channel (rows 1-3 above D, declared: App. A\'s other '
                        'start mode), the gap measured on the inset.  Not reproduced: leg 2\'s side channels beside the inset, '
                        'about 18 um per side wider on the image than the printed single-f schedule draws them (open)'),
               columns={})
    ok_all = True
    for key in ORDER:
        d = DEV[key]
        reg = json.load(open(d['registration']))
        mask = T.Mask(d['dxf'])
        nb = int(reg['measured'].get('neighbour_leg', 2))
        px_um = float(reg['measured']['um_per_px'])             # one image pixel at the inset's 300 um bar
        prof = GA.fig4_profiles(T, raster.im, reg, mask, neighbour=nb)
        if not prof or nb not in prof or 4 not in prof:
            K.fail('%s: the profiles across legs 4 and %d found no walls' % (key, nb))
        if 'metric_search_px' not in reg:
            K.fail('%s: the trace record carries no metric_search_px' % key)
        win_px = float(reg['metric_search_px'])
        nw = GA.narrow_window(T, raster, mask, reg, d['post_max'], leg4_claim(reg), windows=(win_px,))
        walls = reg['measured']['leg2_side_walls']['walls']
        limited = key == 'aligned'
        checks = []

        def check(what, value, limit, claimed=True):
            lim = limited and claimed                 # the limits hold on the claim (leg 4) only
            ok = abs(value) <= limit
            checks.append(dict(what=what, value=value, limit=limit if lim else None, ok=ok if lim else None,
                               claimed=claimed))
            K.say('measure (%s): %-60s %+8.2f um %s' % (key, what, value,
                                                        ('(limit %.2f): %s' % (limit, 'ok' if ok else 'beyond')) if lim
                                                        else ('(reported%s)' % ('' if claimed or not limited else
                                                                                '; not reproduced, register GIF4-1'
                                                                                if abs(value) > limit else ''))))
        for leg in (nb, 4):
            for i, side in enumerate(('low y', 'high y')):
                if prof[leg]['sides_mask'] is None:
                    K.fail('%s: leg %d\'s drawn walls not found on the profiles' % (key, leg))
                check('leg %d side channel (%s), image - drawn' % (leg, side),
                      prof[leg]['sides_image'][i] - prof[leg]['sides_mask'][i], px_um, claimed=(leg == 4))
        for wname, w in walls.items():
            if w['n'] < 20:
                K.say('measure (%s): leg %d wall %s: mean nearest-outline offset %+.2f um over %d samples (fewer than 20: '
                      'listed, not checked)' % (key, nb, wname, w['mean_signed_um'], w['n']))
            else:
                check('leg %d wall %s, mean nearest-outline offset' % (nb, wname), w['mean_signed_um'], px_um,
                      claimed=False)
        nar = nw['%.1f' % win_px]
        check('leg 4 posts\' flanks along the flow, median (the trace\'s window)', nar['post_flanks_along']['median_um'], 0.5)
        check('the claim (leg 4), median |offset| (the trace\'s window)', reg['metric']['median_abs_um'], 1.0)
        rec = json.load(open(d['record']))
        res['columns'][key] = dict(
            mask=d['name'], role=dict(aligned='the package\'s device (the nested reading)', design='the design-level reference',
                                      pathorder='the rejected reading B (shown, not claimed)',
                                      start51='the rejected reading C (shown, not claimed)',
                                      fold20='the 20 um folded-head variant (shown, not claimed)')[key],
            neighbour_leg=nb, image_um_per_px=px_um, profiles=prof, narrow_window=nw, neighbour_walls=walls,
            **{('checks' if limited else 'reference_residuals'): checks},
            trace_metric=dict((k, reg['metric'][k]) for k in ('median_abs_um', 'p90_abs_um', 'n')),
            trace_window_um=reg.get('metric_window_um'), trace_window_px=win_px, window_rule=reg.get('window_rule'),
            null_ratio=K.null_ratio(reg), physics=K.physics_block(rec),
            **({'claim_within': all(c['ok'] for c in checks if c['ok'] is not None)} if limited else {}))
        if limited:
            ok_all = all(c['ok'] for c in checks if c['ok'] is not None)
            res['not_reproduced'] = [c for c in checks if not c['claimed'] and abs(c['value']) > px_um]
    pcs = K.yaml_load(open(DEV['aligned']['preset']))['reproduction']['pc_last_short']
    res['limits_note'] = ('aligned device, its claim (leg 4) only: one image pixel (%.3f um, the inset\'s bar) for leg 4\'s side '
                          'channels; 0.5 um for leg 4\'s posts\' flanks\' median and 1 um for the claim\'s median, at the '
                          'trace\'s own search window (at most half the gap): regression limits, set with the results in view.  '
                          'The neighbouring leg (leg 2) is reported, not limited: its walls and widths are the open item of '
                          'this reproduction (register GIF4-1; not_reproduced lists them).  Not a blind test: the gap (17.2 um) '
                          'was read on this inset\'s posts\' flanks, so the flanks are consistency, and leg 4\'s early end '
                          '(pc_last_short %g um) was read on the inset; the blind quantities are the scale (the bar), the '
                          'rotation (the image), V1 (two independent registrations agree) and the taper direction (V13, '
                          'out/taper_direction.json).  The design, the rejected and the folded-head columns are reported, with '
                          'no limits.' % (px_um, pcs))
    with open(os.path.join(OUT, 'alignment.json'), 'w') as fh:
        json.dump(res, fh, indent=1, default=K.jsonable)
    if not ok_all:
        K.fail('the aligned photomask\'s claim (leg 4) is not within the alignment limits')
    K.say('measure: the aligned photomask\'s claim (leg 4) lies on the inset within its regression limits (set with the '
          'results in view); the neighbouring leg\'s side channels are not reproduced (%d checks beyond one pixel, '
          'register GIF4-1); every column -> out/alignment.json' % len(res.get('not_reproduced', [])))


def _panels(key, T):
    """Return the photomask and its inset panel, after checking that the claim's coloring of the edges is the trace's
    and that the panel was registered on this mask."""
    d = _need_dxf(key)
    if not os.path.exists(d['registration']):
        K.fail('run the trace step first')
    rec = json.load(open(d['registration']))
    mask = T.Mask(d['dxf'])
    claim = leg4_claim(rec)
    x0, x1, y0, y1 = rec['region_um']
    es = mask.edges(x0, x1, y0, y1, post_max=d['post_max'])
    n_claim = sum(1 for e in es if claim(e))
    if (n_claim, len(es) - n_claim) != (rec['claim_edges'], rec['other_edges']):
        K.fail('%s: the claim colouring (%d / %d edges) is not the trace\'s (%d / %d)' % (
            key, n_claim, len(es) - n_claim, rec['claim_edges'], rec['other_edges']))
    K.say('render (%s): leg 4 (the claim, cyan) and leg %s (amber): %d / %d edges of the traced region, as the trace '
          'split them' % (key, rec['measured'].get('neighbour_leg'), n_claim, len(es) - n_claim))
    rasters = {K.sha256(IMAGE): K.load_rgb(IMAGE)}
    panels = K.prepare_panels(PANELS, {'inset': rec}, rasters, mask, claim={'inset': claim})
    for p in panels:
        if p['rec']['mask_sha256'] != mask.sha256:
            K.fail('%s was registered on another mask' % p['key'])
    return mask, panels


def step_render():
    K.say('== render')
    K.use_tree()
    import trace_lib as T
    mask, panels = _panels('aligned', T)
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
    mask_d, panels_d = _panels('design', T)
    zd = K.render_zoom(panels_d, mask_d, os.path.join(DEV['design']['out'], 'fig_zoom_overlay'), overlay=True)
    K.say('render: the reference (design) zoomed panel with its overlay, %d x %d px' % (zd['width_px'], zd['height_px']))
    mask_b, panels_b = _panels('pathorder', T)
    zb = K.render_zoom(panels_b, mask_b, os.path.join(DEV['pathorder']['out'], 'fig_zoom_overlay'), overlay=True)
    K.say('render: the rejected reading B (path order) zoomed panel with its overlay, %d x %d px' % (zb['width_px'], zb['height_px']))
    mask_c, panels_c = _panels('start51', T)
    zc = K.render_zoom(panels_c, mask_c, os.path.join(DEV['start51']['out'], 'fig_zoom_overlay'), overlay=True)
    K.say('render: the rejected reading C (51 um start) zoomed panel with its overlay, %d x %d px' % (zc['width_px'], zc['height_px']))
    mask_f, panels_f = _panels('fold20', T)
    zf = K.render_zoom(panels_f, mask_f, os.path.join(DEV['fold20']['out'], 'fig_zoom_overlay'), overlay=True)
    K.say('render: the 20 um folded-head variant zoomed panel with its overlay, %d x %d px' % (zf['width_px'], zf['height_px']))
    import subprocess                                   # for the V13 figure, which has widths only and no paper raster
    tj = os.path.join(OUT, 'taper_direction.json')
    if not os.path.exists(tj):
        K.fail('run the taper step first')
    items = []
    for key, leg, title in (('pathorder', 'neighbour', 'reading B (path order): the neighbour, leg 3'),
                            ('aligned', 'neighbour', 'the nested reading (the device): the neighbour, leg 2'),
                            ('aligned', 'leg4', 'the control: leg 4\'s array part')):
        items.append('%s=%s' % (os.path.join(DEV[key]['out'], 'taper_%s.json' % leg), title))
    p_ = subprocess.run([sys.executable, '-I', os.path.join(K.TREE, 'tools', 'dissertation_figures.py'), 'v13',
                         os.path.join(OUT, 'fig_taper_direction')] + items, capture_output=True, text=True)
    if p_.returncode:
        K.fail('the V13 figure: %s' % p_.stderr[-300:])
    K.say('render: the V13 figure (reading B, the nested reading, the leg-4 control; widths only) -> out/fig_taper_direction')
    summary = dict(device=DEV['aligned']['name'], paper='Gifford et al. 2014, Lab Chip 14, 4496-4505, Fig. 4 inset',
                   code_stamp=K.STAMP,
                   photomask=dict(file='photomask/%s.dxf' % DEV['aligned']['name'], sha256=mask.sha256),
                   zoom_overlay=z1, zoom_micrographs=z0, full_device=full, full_device_stacked=stack,
                   panels=K.registration_summary(panels),
                   reference_design=dict(photomask=dict(file='design/photomask/%s.dxf' % DEV['design']['name'],
                                                        sha256=mask_d.sha256),
                                         zoom_overlay=zd, panels=K.registration_summary(panels_d)),
                   rejected_pathorder=dict(photomask=dict(file='rejected_pathorder/photomask/%s.dxf' % DEV['pathorder']['name'],
                                                          sha256=mask_b.sha256),
                                           zoom_overlay=zb, panels=K.registration_summary(panels_b),
                                           note='reading B, rejected by the inset\'s taper direction (V13), the caption, the '
                                                'schematic and the arrows (register GIF4-1, INV-2): shown, not claimed'),
                   rejected_start51=dict(photomask=dict(file='rejected_start51/photomask/%s.dxf' % DEV['start51']['name'],
                                                        sha256=mask_c.sha256),
                                         zoom_overlay=zc, panels=K.registration_summary(panels_c),
                                         note='reading C, rejected by the inset\'s beads (register INV-2): shown, not claimed'),
                   variant_fold20=dict(photomask=dict(file='variant_fold20/photomask/%s.dxf' % DEV['fold20']['name'],
                                                      sha256=mask_f.sha256),
                                       zoom_overlay=zf, panels=K.registration_summary(panels_f),
                                       note='the 20 um folded-head variant (spec v5 FOLD start at a 20 um floor): shown, '
                                            'not claimed; not Gifford\'s device'),
                   colours=dict(cyan='the mask\'s edges the panel claims (leg 4)',
                                amber='the mask\'s edges it does not claim (the leg beside leg 4: on the device, leg 2, '
                                      'whose side channels are not reproduced)'),
                   taper_direction='out/taper_direction.json; figure out/fig_taper_direction.{pdf,png}',
                   coordinates='mask um: x along the device, y across; image px: (column, row) of the raster; '
                               'image px = origin + x ex + y ey (map_chamfer)')
    with open(os.path.join(OUT, 'registration.json'), 'w') as fh:
        json.dump(summary, fh, indent=1)


STEPS = dict(check=step_check, build=step_build, trace=step_trace, taper=step_taper, measure=step_measure,
             render=step_render)


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
