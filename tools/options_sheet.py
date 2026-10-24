"""options_sheet.py: generate docs/OPTIONS.md, which gives every option of the generator with its code (YAML key and
generate() argument), its allowed values, its default, its meaning, its glossary symbol, its origin, and the record
field that verifies it on the drawn mask.

The types, allowed values, and defaults are read from opencif/spec.py (SCHEMA, CHOICES, DEFAULTS, DERIVED), so that
the sheet cannot drift from the code, and the descriptions are below. tests/test_options_sheet.py fails if a key or a
value of the schema is undocumented here, or if docs/OPTIONS.md is stale.

Run: python3 tools/options_sheet.py"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from opencif import spec as SP   # noqa: E402
from opencif import api          # noqa: E402

OUT = os.path.join(ROOT, 'docs', 'OPTIONS.md')

# section -> (generate() argument, what the section chooses)
SECTIONS = {
    'name': ('name', 'the device name, which names the files in runs/<name>/'),
    'description': ('description', 'free text, copied to the record'),
    'targets': ('cutoff_um, filtrate_fraction, filtrate_to_retentate_ratio', 'the two design targets (the only '
                                                                             'targets, whereas the depth, gap, posts, '
                                                                             'flow, and fluid are model inputs)'),
    'process': ('process', 'the fabrication process (depth, gap, posts, feature floor, and aspect-ratio limit)'),
    'chip': ('chip', 'the chip area that the layout is fitted to (reported if exceeded and never a pass condition, '
                     'App. B C8 CQ5)'),
    'centre_channel': ('taper, width_in_um, width_min_um', 'the center channel\'s width profile along the device'),
    'layout': ('layout, legs', 'how the device is folded on the chip'),
    'meander': ('meander, seal', 'the meander head that seals the start of the side channels'),
    'inlet': ('inlet', 'the inlet port and what stands between it and gap 1'),
    'outlet': ('outlet', 'the outlet ports and leads'),
    'operation': ('flow_uL_min', 'the operating point'),
    'fluid': ('overrides={"fluid": ...}', 'the fluid (water by default), which enters the Reynolds and Stokes gates '
                                          'only (added in rc7)'),
    'start': ('start', 'the start of the side channels when there is no meander head'),
    'physics': ('overrides={"physics": ...}', 'model options, whose defaults are the validated ones'),
    'gates': ('mode', 'what happens outside the validity envelope'),
    'checks': ('drc, overrides={"checks": ...}', 'validation tolerances'),
    'array': ('array', 'a parallel array of the validated element'),
    'reproduction': ('overrides={"reproduction": ...}', 'engine-level overrides that draw a published device as published'),
    'paper': ('overrides={"paper": ...}', 'free-form source notes, which the engine does not use. Since rc7, '
                                          '`paper.source_trace`, a pointer to the trace of the device against its '
                                          'source figures, is copied to the record (assessment.status.source_trace is '
                                          'NAMED, or NOT NAMED without it), and `paper.targets_fitted`, which names '
                                          'the targets that a reproduction fitted to the figures, is copied to '
                                          'assessment.status.reproduction'),
}

# (section, key) -> (meaning, glossary symbol, verified on the drawn mask by)
KEYS = {
    ('targets', 'cutoff_um'): ('the size cutoff D, above which rigid spheres stay in the center channel; it sets the '
                               'lamina delta* = D/kappa_g',
                               'D = kappa_g delta* (App. A R13(a); implementers step 2)',
                               'Q2 "cutoff D at every gap"; validation.realised_cutoff'),
    ('targets', 'filtrate_fraction'): ('the filtrate fraction F, the fraction of fluid sent to the filtrate (in App. '
                                       'B, "yield" = F)',
                                       'F (App. A R10, eq. 67)', 'Q2 "filtrate fraction F"; validation.network.F_as_drawn'),
    ('targets', 'filtrate_to_retentate_ratio'): ('an alternative to F, the flow ratio r = F/(1-F), so that F = '
                                                 'r/(1+r). It is not the glossary\'s FR (Inglis and Herman\'s '
                                                 'filtration ratio)', 'F/(1-F)', 'as F'),
    ('process', 'depth_um'): ('channel depth h (every duct is h deep)', 'h', 'enters every resistance, and is not '
                                                                             'drawn because a mask is 2-D'),
    ('process', 'gap_um'): ('gap width G between posts, along the flow', 'G (glossary: link width)', 'measured gap openings (Q1)'),
    ('process', 'post_length_um'): ('post size along the flow, which is the wall segment between gaps; the row pitch '
                                    'is Lambda = G + this',
                                    'g_p (glossary: "post width")', 'drawn, and the row pitch enters Q1'),
    ('process', 'post_width_um'): ('post size across the flow, which '
                                   'is the wall that the gap crosses (the gap length)', 'l_g (glossary: "link length")',
                                   'taken from the design rather than measured (see the record notes)'),
    ('process', 'wall_along_flow_um'): ('the glossary\'s name for `post_length_um` (g_p, "post width", the post along '
                                        'the flow), accepted as an alias since rc7; given with `post_length_um`, the '
                                        'two must agree',
                                        'g_p', 'as post_length_um (the record keeps post_length_um)'),
    ('process', 'wall_cross_flow_um'): ('the glossary\'s name for `post_width_um` (l_g, "link length", the wall that '
                                        'the gap crosses), accepted as an alias since rc7; given with `post_width_um`, '
                                        'the two must agree',
                                        'l_g', 'as post_width_um (the record keeps post_width_um)'),
    ('process', 'min_feature_um'): ('the '
                                    'process floor, which is the narrowest channel or wall that can be made', 'F_min',
                                    '"DRC at the process floor", "side channels continuous", "no intersections"'),
    ('process', 'min_radius_um'): ('the smallest inner-corner radius that the '
                                   'process makes', 'R_min (glossary, Validity-Envelope and External Inputs: external)',
                                   '"corner radius"'),
    ('process', 'fillet_um'): ('fillet radius (closing, then opening), drawn over only two zones, the inlet with '
                               'the meander head, and the outlet (its lead steps and port necks). The linear legs are not '
                               'filleted, because their posts are stadiums and the U-turns and hairpins are arcs. 0 '
                               'draws square corners, as the papers\' CAD drawings do, whereas Dinh\'s fabricated head '
                               '(Fig. 2B) shows its corners rounded (about 4-8 um on its 50 um bar)',
                               'r_f', 'drawn; enters the head T-junction term (S_T) and its band; "corner radius"'),
    ('process', 'gap_entrance_c'): ('gap entrance coefficient in R_gap = r(G) [l_g + c_ent D_h]', 'c_ent (glossary, Validity-Envelope and External Inputs: external)',
                                    'a design value (future work: 3-D calibration)'),
    ('process', 'aspect_ratio_max'): ('the deepest channel that the process makes, as depth over width (AR_max). '
                                      'CIF_v2\'s floor w_min = max(min_feature, depth/AR_max) is a gate, checked '
                                      'before drawing, on the gap, the meander channel, every center width, and every '
                                      'straight side channel (spec v5 G7; G11 would fold a narrower side channel, '
                                      'which opencif does not draw), and on the drawn mask (the narrowest CIF '
                                      'channel\'s depth/width). Strict mode refuses or fails, whereas evaluate and '
                                      'reproduce modes report (LIMITATIONS 7). The default is CIF_v2\'s 10 (Gifford '
                                      '2014 p. 4503 assumes "a 10 : 1 aspect ratio limitation"), and a finer process '
                                      'is a stated capability rather than a foundry fact',
                                      'AR_max (glossary, Validity-Envelope and External Inputs: external); w_min (spec v5 sec. 18)',
                                      '"aspect ratio of the narrowest CIF channel" (a pass condition in strict mode, reported in '
                                      'evaluate and reproduce modes); assessment.over_aspect_ratio'),
    ('chip', 'width_mm'): ('chip width', '', '"fits the chip area" (reported rather than a pass condition; '
                                             'assessment.over_chip_area and Result.over_chip_area say so, and '
                                             'Result.ok stays true); sets the head-width search range'),
    ('chip', 'length_mm'): ('chip length', '', '"fits the chip area" (reported rather than a pass condition); '
                                               'leg-count fitting'),
    ('chip', 'margin_um'): ('keep-out margin at the chip edge', '', '"fits the chip area" (reported rather than a pass '
                                                                    'condition)'),
    ('centre_channel', 'taper'): ('the center-width profile (values below)', 'w_c(i)', 'measured center widths (Q1)'),
    ('centre_channel', 'width_in_um'): ('center width '
                                        'at gap 1', 'w_c(0)', 'measured (validation.drawn_parameters.wc_in_meas)'),
    ('centre_channel', 'width_min_um'): ('narrowest '
                                         'center width (constant_shear only)', 'w_c,min', 'measured (wc_out_meas)'),
    ('centre_channel', 'shear_formula'): ('constant_shear only, the wall shear that the taper holds constant (values '
                                          'below; added in rc7)',
                                          'gamma_w', 'measured center widths (Q1); the record\'s design.wc'),
    ('layout', 'style'): ('how the legs are arranged (values below)', '', 'components_drawn.legs, u_turns'),
    ('layout', 'legs'): ('number of legs, or "auto" for the fewest that fit the chip with no intersections', '',
                         'components_drawn.legs; "drawn components as requested"'),
    ('layout', 'order'): ('paperclip only, the legs nested (the default) or side by side in path order (values below; '
                          'added in rc7)',
                          '', 'design.leg_Y; "no intersections"'),
    ('layout', 'leg_spacing_um'): ('the U-turn radius for a serpentine and the leg spacing for a paperclip. The '
                                   'default is max(2 x the widest row\'s half-width, the head\'s half-width + it) + '
                                   '150 um for a paperclip and that half-width + floor/2 + 20 um for a serpentine, '
                                   'each grown until the loop channels and the disc inside them are above the floor. A '
                                   'given value is used as given, and parts that then collide are refused before '
                                   'anything is written (rc7)', 'R_turn / d', '"no intersections"'),
    ('meander', 'enabled'): ('draw a meander head (true) or start the side channels at gap 1 (false). The YAML default '
                             'is false, and the generate() default is a packing head (meander=\'pack\'; meander=None '
                             'for none)', '',
                             'components_drawn.head_units'),
    ('meander', 'style'): ('the head style (values below)', '', 'components_drawn.head_drawn, head_arms, arm_lane_pairs, '
                                                               'head_unit_widths_um'),
    ('meander', 'arms'): ('arms style only, the number of arms (Dinh: 3)', '', 'components_drawn.head_arms'),
    ('meander', 'channel_um'): ('width of the head\'s units (its channel), GS, where GS/G must lie in gamma_GS [1.4, '
                                '1.5] (spec v5 sec. 14.4, 18; strict mode refuses a value outside it, and the other '
                                'modes report it)', 'G_S',
                                'components_drawn.head_unit_widths_um; gamma_GS pre-flight'),
    ('meander', 'wall_um'): ('wall between neighboring head units', '', 'DRC'),
    ('meander', 'narrow_um'): ('narrow style only, the width of the narrow units', 'w_narrow', 'head_unit_widths_um'),
    ('meander', 'head_width_um'): ('lateral room for the head, or "auto" (the narrowest that seals gap 1)', 'W',
                                   '"head inside its envelope"'),
    ('meander', 'seal'): ('whether gap 1 must extract exactly f* (values below)', 'f(1) = f*', 'components_drawn.gap1_sealed, '
                                                                                         'gap1_f_over_fstar'),
    ('meander', 'transition'): ('how the last head unit hands over to the straight side channel (values below)', '',
                                'Q1 at the seam rows'),
    ('meander', 'rectilinear_rows'): ('wall-block rows between the head and the first stadium post (>= 1)', 'n_rect', 'topology'),
    ('meander', 'seam_side_width_um'): ('side-channel width required at the seam', 'w_s(in)', 'Q1 at the seam'),
    ('meander', 'seam_at'): ('where seam_side_width_um is required (values below)', '', 'Q1 at the seam'),
    ('meander', 'turns'): ('how the head\'s quarter turns are drawn (values below; added in rc7), whereas the U-turns '
                           'stay arcs', '',
                           'Q1 on the head units; "corner radius"'),
    ('meander', 'corner_radius_um'): ('square turns only, the rounding of each square corner\'s outer (convex) corner, '
                                      'from 0 to 0.307 of the narrowest head unit (the 2-D solves\' range, refused '
                                      'beyond it), drawn and modeled at no less than the head fillet, because the head '
                                      'smoothing rounds every corner at it. Within that range it has no measurable '
                                      'effect on the corner term (added in rc7)',
                                      'r_out', 'drawn; "corner radius"'),
    ('meander', 'junction_radius_um'): ('the rounding of the wall corners where each head unit\'s stem leaves the '
                                        'center channel (the head gaps\' entrance), in [0, half the thinner of '
                                        'post_width_um and wall_um], for the drawing only, since the gap entrance term '
                                        'c_ent is an external input (added in rc7)', 'r_j',
                                        'drawn ("corner radius"); the head gaps are measured mid-wall, above the roundings'),
    ('meander', 'inner_corner_radius_um'): ('square turns only, the fillet at each square corner\'s inner corner (the '
                                            'end of the wall between lanes), from 0 to half the narrower of the '
                                            'narrowest head unit and wall_um, drawn and modeled at no less than the '
                                            'head fillet. It enters the corner term (schedule.c_corner, interpolated '
                                            'in r_in/w; added in rc7)', 'r_in',
                                            'drawn; Q1 (the measurement takes the roundings\' area off, because the '
                                            'corner term counts it)'),
    ('inlet', 'style'): ('what stands between the inlet port and gap 1 (values below)', '', 'components_drawn.filter_posts_counted'),
    ('inlet', 'port_radius_um'): ('inlet port radius', '', '"inlet and outlet ports apart"'),
    ('inlet', 'stem_half_width_um'): ('half width of the inlet stem', '', 'DRC, topology'),
    ('inlet', 'port_gap_um'): ('gap from the inlet port\'s edge to the head\'s top (default 150 um)', '', '"no intersections"'),
    ('outlet', 'style'): ('outlet arrangement (values below)', '', 'components_drawn.outlet_style, outlet_ports'),
    ('outlet', 'port_radius_um'): ('filtrate port radius', '', '"outlet ports separated"'),
    ('outlet', 'centre_port_radius_um'): ('retentate (center) port radius', '', '"outlet ports separated"'),
    ('outlet', 'stretch_um'): ('straight outlet stretch after the last gap. A serpentine or paperclip whose last leg '
                               'ends beside a loop replaces it by the stretch that clears the loop, so the spec\'s '
                               'value is not used there, and the key is refused with `flare_deg`, which sets the '
                               'stretch. Since rc5, where the side channel meets a lead of another width (every '
                               'teardrop and ring outlet, and a separate outlet only when w_s(N) < 40 um), the '
                               'change\'s excess over the 1-D rule enters the closure as a junction term '
                               '(schedule.outlet_exp_c, from 2-D Stokes and Hele-Shaw solves, depth-weighted and '
                               'banded)', 'c_end',
                               'validation.outlet_measured (stretch_end); the model band (junction band)'),
    ('outlet', 'flare_deg'): ('a flared outlet lead (added in rc5), with the side channel\'s outer wall diverging at '
                              'this angle per side (1-45) from w_s(N) to the lead width after `flare_start_um`, and '
                              'its resistance integrated (exact wall rule) plus the flare\'s excess over that integral '
                              '(schedule.outlet_exp_c, as for a step). It applies to teardrop and separate outlets and '
                              'is refused with the ring, with `stretch_um`, and with `physics.fail_safe_margin`. A '
                              'lead no wider than w_s(N) draws no flare and keeps the stretch that it has without one, '
                              'and an unset value gives the right-angle step at the end of `stretch_um`',
                              '', 'validation.outlet_measured (the stretch integrated); the model band'),
    ('outlet', 'flare_start_um'): ('flare_deg only, the straight run from the last gap\'s center to the flare\'s '
                                   'onset, which is the virtual corner of its 25 um arc (Dinh 2024 Fig. 2B, fourth '
                                   'panel: 65.4 +- 2 um, measured at the 50 um bar read with its anti-aliased ends; '
                                   'rc5 used 67). Where the last leg ends beside a loop, the run is the longer of this '
                                   'and the loop\'s clearance',
                                   '', 'validation.outlet_measured; design.outlet (L_out_um, stretch_end)'),
    ('outlet', 'tubing'): ('off-chip tubing on the outlets (centre_id_mm, side_id_mm, length_m), which enters the '
                           'outlet closure',
                           'R_ext', 'reported in the tubing scenarios of the model band'),
    ('outlet', 'ring'): ('ring geometry (radius_um, centre_offset_um, filtrate_port_offset_um, turn_radius_um), scaled '
                         'on the ports by default', '', 'outlet_measured'),
    ('operation', 'flow_uL_min'): ('total inlet flow (the fluid is set in section fluid, water by default)', 'Q_in',
                                   'the validity envelope, with channel Re <= 10 and the Stokes gates Re_gap and '
                                   'Re_delta <= 1e-2'),
    ('start', 'style'): ('no-head start rule (values below)', 'sigma_0', 'Q2 (a floor start\'s entrance rows are declared in reproduce mode; outside it they usually FAIL Q2)'),
    ('start', 'side_width_um'): ('floor start '
                                 'only, the first side-channel width (the floor by default)', 'w_s(1)', 'measured'),
    ('physics', 'fail_safe_margin'): ('lower f on every row by the largest over-extraction that any bias scenario '
                                      'gives '
                                      '(Xia\'s pull-in made quantitative); off by default', '', 'robustness scenarios'),
    ('physics', 'trim_F'): ('trim the last gap so that F = target exactly; off by default', '', 'Q2 F'),
    ('physics', 'wall_rule'): ('the straight side-wall rule, exact (resistance-exact, the default) or legacy (the old '
                               'drawing, for comparison)',
                               '', 'Q1'),
    ('physics', 'bend_model'): ('stokes for the exact curved-duct correction at every arc (the default), or none for '
                                'the plain 1-D rule, for comparison only',
                                '', 'model band'),
    ('gates', 'mode'): ('what happens outside the validity envelope (values below)', '', 'assessment.mode, verdict, exit code'),
    ('checks', 'drc'): ('run the DRC '
                        '(true by default), without which the verdict is INCOMPLETE', '', '"DRC at the process floor"'),
    ('checks', 'cutoff_tolerance'): ('Q2 tolerance on every gap\'s cutoff, as a fraction of D (default 0.01, at most 0.05)', '',
                                     'Q2'),
    ('checks', 'F_tolerance'): ('Q2 tolerance on F (default 1e-3, at most 0.01)', '', 'Q2'),
    ('checks', 'passage_report'): ('report each measured channel class (gap openings, straight side channels, U-turn '
                                   'sides, head units, and outlet leads, with the center channel as reference) against '
                                   'D and 2 D, with a warning note below 2 D, and, where rows are realized above D, '
                                   'against the largest particle that a gap passes, min(D_realised_max, G), with its '
                                   'own warning (added in rc7). It is a report rather than a check, because opencif '
                                   'has no passage rule (LIMITATIONS 26)', '', 'validation.passage'),
    ('fluid', 'density_kg_m3'): ('the fluid\'s '
                                 'density (added in rc7)', 'rho', 'the Reynolds gates (validation.gates.fluid)'),
    ('fluid', 'viscosity_Pa_s'): ('the fluid\'s viscosity (added in rc7). The schedule does not depend on it, because '
                                  'the resistance ratios are independent of mu for a Newtonian fluid, whereas the '
                                  'gates do', 'mu',
                                  'the Reynolds and Stokes gates (validation.gates.fluid)'),
    ('fluid', 'rheology'): ('the fluid model (values below; '
                            'added in rc7), of which only a Newtonian fluid is modeled', '', 'validation.gates.fluid'),
    ('array', 'count'): ('number of copies of the validated element', 'n', 'array record'),
    ('array', 'pitch_um'): ('spacing '
                            'between copies (by default the element\'s height plus clearance)', '', 'array DRC'),
    ('array', 'inlet'): ('shared bus or separate inlets (values below)', '', 'array record'),
    ('array', 'bus_width_um'): ('width of the shared inlet bus', '', 'array record (flow split, Reynolds at the busiest element)'),
}

# (section, key, value) -> (what is drawn, origin, example spec, restrictions)
VALUES = {
    ('centre_channel', 'taper', 'constant'): ('w_c constant; f constant along the device', 'Gifford 2014',
                                              'examples/gifford_physics_linear.yaml', 'width_min_um not allowed'),
    ('centre_channel', 'taper', 'constant_shear'): ('w_c narrows from width_in_um to width_min_um as the flow drops, '
                                                    'holding the wall shear, and '
                                                    'f rises to hold delta*', 'Xia 2016 (Dinh 2024\'s form assumed)',
                                                    'examples/matrix/m02_linear_shear_arms_filter_teardrop.yaml',
                                                    'width_min_um required'),
    ('layout', 'style', 'linear'): ('one straight leg', 'Gifford 2014 Fig. 1', 'examples/matrix/m01_linear_const_pack_plain_separate.yaml',
                                    'legs = 1'),
    ('layout', 'style', 'serpentine'): ('legs side by side joined by U-turns; the outlet leaves away from the head',
                                        'this project', 'examples/matrix/m07_serpentine_const_pack_filter_separate.yaml',
                                        'legs odd (>= 3); no ring outlet'),
    ('layout', 'style', 'paperclip'): ('nested legs (each loop wraps the previous)', 'Gifford 2014 Fig. 4, Xia 2016 Fig. 3a',
                                       'examples/matrix/m11_paperclip_const_pack_filter_ring.yaml', 'legs >= 2'),
    ('layout', 'order', 'nested'): ('DEFAULT. Leg j at Y = 0, +d, -d, +2d, ..., so each loop wraps the previous one, '
                                    'all turning '
                                    'alike', 'Gifford 2014 Fig. 4 schematic, caption and inset (RC8), Xia 2016 Fig. 3a',
                                    'presets/gifford2014_fig4_aligned.yaml', ''),
    ('layout', 'order', 'path'): ('Added in rc7. Leg j at Y = j d (or a reproduction.leg_Y stepping one way), so the '
                                  'legs lie side by side in path order and the loops alternate their turning '
                                  'direction; the outlet leaves at the inlet end for an even leg '
                                  'count', 'Gifford 2014 Fig. 4 inset, reading B (register GIF4-1; rejected in RC8)',
                                  'provenance/rc7_traces/variants/gifford2014_fig4_pathorder.yaml', 'paperclip only'),
    ('meander', 'style', 'pack'): ('DEFAULT. Dinh\'s arms topology with as many arms as the head width holds, folded '
                                   'into lanes, then hairpins, '
                                   'at the narrowest width that seals gap 1', 'this project (after Dinh 2024 and the '
                                   'student generator\'s packing)', 'examples/matrix/m11_paperclip_const_pack_filter_ring.yaml', ''),
    ('meander', 'style', 'arms'): ('Dinh\'s head, with `arms` single-lane arms and then hairpins', 'Dinh 2024 Fig. 2A',
                                   'examples/matrix/m12_paperclip_shear_arms_plain_ring.yaml', ''),
    ('meander', 'style', 'hairpins'): ('Xia\'s head, '
                                       'with perpendicular hairpin units only, all of one width', 'Xia 2016 Fig. 1c',
                                       'examples/matrix/m03_linear_const_hairpins_filter_ring.yaml',
                                       'cannot seal gap 1 at most targets, so it needs seal: allow (evaluate mode), '
                                       'and gap 1 is declared'),
    ('meander', 'style', 'narrow'): ('opt-in, with process-floor units first and then full-width units (a smaller head '
                                     'with about three times the width '
                                     'sensitivity)', 'this project', 'examples/matrix/m09_paperclip_const_narrow_plain_teardrop.yaml',
                                     ''),
    ('meander', 'seal', 'require'): ('DEFAULT. Gap 1 at f* (its cutoff is D), or the build is REFUSED with the numbers '
                                     '(the arms needed against the room, and the width that would seal), and no other '
                                     'head is substituted', '', '', ''),
    ('meander', 'seal', 'auto'): ('as require in strict mode; otherwise an unsealed gap 1 is allowed and declared', '', '', ''),
    ('meander', 'seal', 'allow'): ('the head asked for may leave gap 1 unsealed (extracting more, with a cutoff above '
                                   'D), and gap 1 is declared',
                                   '', 'examples/matrix/m10_paperclip_shear_hairpins_filter_separate.yaml', 'not in strict mode'),
    ('meander', 'transition', 'exact'): ('DEFAULT. Rectilinear rows at their solved widths', '', '', ''),
    ('meander', 'transition', 'paper'): ('the papers\' transition, in which the side channel beside the wall blocks is '
                                         'G_S wide and narrows at the first '
                                         'post (their deliberate pull-in)', 'Dinh 2024 Fig. 2B, Xia 2016 Fig. 1c(ii)',
                                         'presets/dinh2024_v1_specified.yaml', 'reproduce mode only'),
    ('meander', 'seam_at', 'segment'): ('DEFAULT. '
                                        'seam_side_width_um at the first straight segment after the head', '', '', ''),
    ('meander', 'seam_at', 'post'): ('seam_side_width_um at the first stadium post, after the rectilinear rows',
                                     'Xia 2016 Fig. 1c(ii)', 'presets/xia2016_rs7_specified.yaml', ''),
    ('meander', 'turns', 'arc'): ('DEFAULT. Each quarter turn is an arc of centerline radius s/2 (rc6)', '', '', ''),
    ('meander', 'turns', 'square'): ('Added in rc7. Each quarter turn is a right-angle corner, with its centerline '
                                     'through the corner point (2r of path where the arc has pi r/2), the corner term '
                                     'c_corner w (from 2-D Hele-Shaw and Stokes solves, depth-weighted and banded), '
                                     'and the roundings `corner_radius_um` '
                                     'and `inner_corner_radius_um`', 'Dinh 2024 Fig. 2B (register DINH-1)',
                                     'presets/dinh2024_v1_depicted.yaml', 'a meander head only'),
    ('centre_channel', 'shear_formula', 'exact'): ('DEFAULT. The rectangular duct\'s maximum wall shear C(alpha) Q/(a '
                                                   'b^2) (CIF_v2 C_shear)', '', '', 'constant_shear only'),
    ('centre_channel', 'shear_formula', 'parallel_plates'): ('Added in rc7. The parallel-plate shear 6 Q/(a b^2), '
                                                             'which makes w_c proportional to the center flow. It is a '
                                                             'departure from spec v5 S16-b (constant shear with '
                                                             'C(alpha); "NEVER 6") and is the form that Xia 2016\'s '
                                                             'Fig. 3 panels d and e read, since Xia\'s text does not '
                                                             'state the coefficient. The record reports the true peak '
                                                             'shear through C(alpha) (gates.shear_peak)',
                                                            'Xia 2016 Fig. 3 (register XIA-1)',
                                                            'presets/xia2016_rs7_aligned.yaml', 'constant_shear only; reproduce mode only'),
    ('fluid', 'rheology', 'newtonian'): ('DEFAULT and the only value, a Newtonian fluid', '', '', ''),
    ('inlet', 'style', 'plain'): ('DEFAULT. A plain stem from '
                                  'the port to gap 1', '', 'examples/matrix/m01_linear_const_pack_plain_separate.yaml',
                                  'no posts may be drawn (checked)'),
    ('inlet', 'style', 'filter'): ('an aggregate filter, which is a bulb (half-width up to 150 um) of staggered 30 um '
                                   'posts at 60 um pitch, between the arms or above the head, upstream of gap 1 and '
                                   'outside the network', 'this project',
                                   'examples/matrix/m11_paperclip_const_pack_filter_ring.yaml', 'posts counted on the drawing'),
    ('inlet', 'style', 'dinh_filter'): ('Dinh\'s two-stage filter, in its own absolute coordinates', 'Dinh 2024 Fig. 2A',
                                        'presets/dinh2024_v1_specified.yaml', 'only with meander style arms, arms 3'),
    ('outlet', 'style', 'separate'): ('DEFAULT. Each side channel goes to its own filtrate port, and the center '
                                      'channel to the retentate port',
                                      'Gifford 2014 Fig. 4', 'examples/matrix/m07_serpentine_const_pack_filter_separate.yaml',
                                      '3 ports'),
    ('outlet', 'style', 'teardrop'): ('both '
                                      'sides wrap around the center port into one filtrate port', 'Dinh 2024 Fig. 2A',
                                      'examples/matrix/m06_serpentine_shear_pack_filter_teardrop.yaml', '2 ports'),
    ('outlet', 'style', 'ring'): ('both side leads follow one circle around the center port into a filtrate port '
                                  'beyond it',
                                  'after Xia 2016 Fig. 3a', 'examples/matrix/m11_paperclip_const_pack_filter_ring.yaml',
                                  'linear and paperclip layouts only; REFUSED (exit 5, any mode) when the ring cannot '
                                  'be closed at any leg count that the layout tries (plan-stage lattice, three of the '
                                  'four processes at F 0.9; tests/test_lattice.py)'),
    ('start', 'style', 'sealed'): ('DEFAULT. The side channel starts sealed (w_s(0) = 0), with gap 1 at f*, which '
                                   'needs a fine process without a '
                                   'head', 'Gifford 2014 Fig. 2', 'examples/matrix/m13_paperclip_const_nomeander_sealed_fine.yaml',
                                   'no-head devices'),
    ('start', 'style', 'floor'): ('the first side channel drawn at the floor (side_width_um), so the entrance '
                                  'over-extracts; the '
                                  'design network predicts those rows, and they are declared', 'Gifford 2014 p. 4498',
                                  'examples/matrix/m14_paperclip_const_nomeander_floorstart19.yaml', 'no-head devices only'),
    ('physics', 'wall_rule', 'exact'): ('DEFAULT. Straight side walls drawn so that their integrated resistance is the '
                                        'closure\'s', '', '', ''),
    ('physics', 'wall_rule', 'legacy'): ('the older '
                                         'linear wall rule (errors up to 15% in f), for comparison only', '', '', ''),
    ('physics', 'bend_model', 'stokes'): ('DEFAULT. Exact '
                                          'Stokes curved-duct correction less its junctions at every arc', '', '', ''),
    ('physics', 'bend_model', 'none'): ('the plain 1-D rule '
                                        'at bends, for comparison only', '', '', 'not in strict mode; not arrayable'),
    ('gates', 'mode', 'strict'): ('DEFAULT. Outside the validity envelope (the kappa plateau and foot, channel Re, and '
                                  'the Stokes gates from the App. B thesis; gap dominance M13, thin lamina E11, and '
                                  'taper rate M16 from CIF_v2), or failing a gate checked before drawing (a meander '
                                  'head whose GS/G is outside gamma_GS [1.4, 1.5], or the gap, a drawn GS, a center '
                                  'width, a straight side channel, or the terminal side width below w_min = max(floor, '
                                  'depth/AR_max), from spec v5 sec. 14.4, 18, G7 and CIF_v2 design()), the build is '
                                  'REFUSED before drawing (exit 5), and the record still gives the cutoff band. On the '
                                  'drawn mask, the narrowest CIF channel above the aspect-ratio limit is a FAIL (exit '
                                  '2) in strict mode and is reported in the other modes (LIMITATIONS 7). A design over '
                                  'its chip area is drawn and reported, and it is never refused or failed (App. B C8 '
                                  'CQ5)',
                                  'App. B thesis; CIF_v2 M13/M14/E11/M16; spec v5 sec. 14.4, 18-19, G7', '', ''),
    ('gates', 'mode', 'evaluate'): ('the mask is drawn and fully validated, and leaving the envelope gives exit 3 '
                                    '(`opencif check` also exits 3 on the plan); a pre-flight gate is reported in the '
                                    'notes', '',
                                    'examples/dinh_targets_sealed_head.yaml', ''),
    ('gates', 'mode', 'reproduce'): ('as evaluate, and the papers\' deliberate departures (reproduction block, transition: paper) '
                                     'are allowed, each declared and attributed row by row', '', 'presets/*.yaml', ''),
    ('array', 'inlet', 'shared'): ('one inlet bus feeds every copy (drawn only when every inlet can be reached in one layer)',
                                   '', 'tools/make_arrays.py m07_... 3 shared 25 75', ''),
    ('array', 'inlet', 'separate'): ('each copy keeps its own inlet port', '', '', ''),
}

# defaults that live in the engine or the pipeline rather than in spec.DEFAULTS (checked by tests/test_options_sheet.py)
EXTRA_DEFAULTS = {('layout', 'legs'): '1 (linear); `auto` otherwise: the fewest that fit the chip, or 4 (paperclip) / 5 '
                                      '(serpentine) without a chip area',
                  ('layout', 'leg_spacing_um'): 'auto',
                  ('meander', 'arms'): '`3`', ('inlet', 'port_gap_um'): '`150`', ('outlet', 'stretch_um'): '`200`',
                  ('outlet', 'flare_deg'): 'unset (the step)', ('outlet', 'flare_start_um'): '`65.4`',
                  ('outlet', 'tubing'): 'none', ('outlet', 'ring'): 'scaled on the ports',
                  ('start', 'side_width_um'): 'the floor (`min_feature_um`)',
                  ('meander', 'junction_radius_um'): '`0`',
                  ('checks', 'drc'): '`True`', ('checks', 'cutoff_tolerance'): '`0.01`', ('checks', 'F_tolerance'): '`0.001`',
                  ('array', 'pitch_um'): 'auto', ('array', 'inlet'): '`shared`', ('array', 'bus_width_um'): 'auto'}

PROCESS_SOURCES = {
    'dinh2024': '`presets/dinh2024_v1_specified.yaml`, with G 19 um and depth 140 um printed (Dinh 2024 p. 917, Fig. 2 '
                'caption) and the post sizes measured on Fig. 2B (PRESET_PROVENANCE.md)',
    'xia2016': '`presets/xia2016_rs7_specified.yaml`, Xia 2016 RS-7, with G 16 um printed, depth 140 um inside the '
               'printed 135-150 um, and the post sizes measured on Fig. 3g (PRESET_PROVENANCE.md)',
    'gifford2014': '`presets/gifford2014_fig4_design.yaml` with a 19 um floor, as `examples/gifford_targets_meander_head_19um.yaml`',
    'generic19': '`examples/matrix/*.yaml`, this project\'s 19 um, 100 um deep process',
}
GEN_ARGS = [
    ('`cutoff_um`', '`targets.cutoff_um`', 'D in um', '**required** (1st argument)'),
    ('`filtrate_fraction`', '`targets.filtrate_fraction`', 'F in (0, 1)', '**required** (2nd argument), or:'),
    ('`filtrate_to_retentate_ratio`', '`targets.filtrate_to_retentate_ratio`', 'F/(1-F) > 0', 'none'),
    ('`process`', '`process`', 'a name below or a dict', '`generic19`'),
    ('`chip`', '`chip`', '(width_mm, length_mm[, margin_um]) or a dict', 'none (no area report and no leg fitting)'),
    ('`layout`, `legs`', '`layout.style`, `layout.legs`', 'linear / serpentine / paperclip; int or `auto`', '`linear`; legs as section 3'),
    ('`taper`', '`centre_channel`', 'constant / constant_shear, or (style, width_in_um[, width_min_um])', '`constant`'),
    ('`width_in_um`', '`centre_channel.width_in_um`', 'um', 'D / %.2f, so u = D/w_c = %.2f, '
                                                            'on the kappa plateau' % (api.DEFAULT_U, api.DEFAULT_U)),
    ('`width_min_um`', '`centre_channel.width_min_um`', 'um (constant_shear only)', '0.7 x width_in_um'),
    ('`meander`', '`meander`', 'pack / arms / hairpins / narrow, a dict, or None (no head)',
     '`pack`, the packing head, which is this project\'s own design (no paper has it, and its seam rows carry the 2-D '
     'junction fits, LIMITATIONS 8); the YAML default is no head'),
    ('`seal`', '`meander.seal`', 'require / auto / allow', '`require`'),
    ('`inlet`', '`inlet`', 'plain / filter / dinh_filter, or a dict', '`plain`'),
    ('`outlet`', '`outlet`', 'separate / teardrop / ring, or a dict', '`separate`'),
    ('`start`', '`start`', 'sealed / floor, or a dict (no head only)', '`sealed`'),
    ('`flow_uL_min`', '`operation.flow_uL_min`', 'uL/min', '`30`'),
    ('`mode`', '`gates.mode`', 'strict / evaluate / reproduce', '`strict`'),
    ('`array`', '`array`', 'an int (count) or a dict', 'none'),
    ('`name`, `description`', '`name`, `description`', 'text', '`cif_D<D>_F<F>` (give `name=` when two designs share D and F)'),
    ('`overrides`', 'any section', 'a dict of section dicts', 'none'),
    ('`out`, `render`, `drc`, `check_only`, `verbose`', '-', 'output folder; figures; DRC; check_only=True for the '
                                                             'plan, the gates, and the outline\'s collision checks in '
                                                             'memory, with nothing written (\'plan\' for the plan and '
                                                             'gates '
                                                             'only); log', '`runs/<name>`, True, True, False, False'),
]

REQUIRED = [('targets', 'cutoff_um'), ('targets', 'filtrate_fraction'), ('process', 'depth_um'), ('process', 'gap_um'),
            ('process', 'min_feature_um'), ('centre_channel', 'width_in_um')]

MATRIX_CODES = [
    ('layout', 'linear / serpentine / paperclip', 'layout.style'),
    ('taper', 'const = constant; shear = constant_shear', 'centre_channel.taper'),
    ('head', 'pack / arms / hairpins / narrow; nomeander = meander.enabled: false', 'meander.style'),
    ('inlet', 'plain / filter', 'inlet.style'),
    ('outlet', 'separate / teardrop / ring', 'outlet.style'),
    ('suffix', 'sealed_fine = start.style sealed on a 4.3 um process; floorstart19 = start.style floor at 19 um', 'start.style'),
]


def _fmt(v):
    """Return a schema type as the sheet prints it, with 'mapping' for a nested section."""
    if isinstance(v, tuple):
        return '/'.join(t.__name__ for t in v)
    if isinstance(v, type):
        return v.__name__
    return 'mapping'


def leaves():
    """Return (section, key, type) for every leaf of the schema, with a dict-valued key documented as one row."""
    out = []
    for sec, sub in SP.SCHEMA.items():
        if isinstance(sub, dict):
            for k, t in sub.items():
                out.append((sec, k, t))
        else:
            out.append((sec, None, sub))
    return out


def default_of(sec, key):
    """Return the default that the sheet states for a key, from the spec's defaults, its derived values, the required
    keys, or EXTRA_DEFAULTS."""
    if (sec, key) in SP.DEFAULTS:
        return '`%s`' % (SP.DEFAULTS[(sec, key)],)
    if (sec, key) in SP.DERIVED:
        return 'derived: %s' % SP.DERIVED[(sec, key)]
    if (sec, key) in REQUIRED:
        return '**required**'
    if (sec, key) in EXTRA_DEFAULTS:
        return EXTRA_DEFAULTS[(sec, key)]
    if (sec, key) == ('centre_channel', 'width_min_um'):
        return 'required for `constant_shear`'
    if (sec, key) == ('targets', 'filtrate_to_retentate_ratio'):
        return 'instead of F'
    return ''


def render():
    """Return the text of docs/OPTIONS.md."""
    L = []
    w = L.append
    w('# opencif options sheet')
    w('')
    w('This sheet is generated by `tools/options_sheet.py` from `opencif/spec.py` (the schema, the allowed values, and '
      'the defaults), and `tests/test_options_sheet.py` fails if an option is undocumented or this file is stale.')
    w('')
    w('## 1. What the generator needs')
    w('')
    w('The two design targets are the size cutoff D and the filtrate fraction F. The other inputs are the fabrication '
      'context (process and chip) and the parts (layout, taper and center width, meander head, inlet, and outlet). The '
      'gap fraction f_gap and the row count N are outputs rather than inputs, except in reproduce mode, where a '
      '`reproduction` block may prescribe a paper\'s depicted head fractions and those rows are declared. On the kappa '
      'plateau, D and F are independent (App. A R13(c); App. B C8).')
    w('')
    w('```python')
    w('import opencif')
    w("r = opencif.generate(4.0, 0.5, out='/tmp/cif')    # D = 4 um, F = 0.5; every other choice defaulted")
    w("r = opencif.generate(2.5, 0.6, process='dinh2024', chip=(25, 75), layout='serpentine',")
    w("                     taper=('constant_shear', 70, 50), meander='pack', inlet='filter', outlet='teardrop')")
    w('print(r.verdict, r.exit_code, r.dxf, r.figures)   # r.record is the full validation record')
    w("opencif.generate(4.0, 0.5, check_only=True)        # plan, gates, outline collisions; nothing written")
    w('```')
    w('')
    w('```')
    w('python3 -m opencif generate 2.5 0.6 --process dinh2024 --chip 25x75 --layout serpentine \\')
    w('        --taper constant_shear --width-in 70 --width-min 50 --meander pack --inlet filter --outlet teardrop')
    w('python3 -m opencif generate 4 0.5 --check                      # plan, gates, outline collisions; nothing written')
    w('python3 -m opencif generate 4 0.5 --set checks.cutoff_tolerance=0.02 --set physics.trim_F=true')
    w('python3 -m opencif build runs/<name>/<name>_spec.yaml           # rebuild from the spec generate() wrote')
    w('```')
    w('')
    w('### generate() arguments and their defaults')
    w('')
    w('| argument | YAML section / key | accepts | default in generate() |')
    w('|---|---|---|---|')
    for row in GEN_ARGS:
        w('| %s | %s | %s | %s |' % row)
    w('')
    w('Every part argument also takes a dict of its section\'s keys (section 3), for example `meander={\'style\': '
      '\'arms\', \'arms\': 3, \'seal\': \'require\'}`. A `process=` dict replaces the named process entirely, so unset '
      'post sizes then take the spec defaults, 25 and 50 um. `overrides=` merges a dict into any section last, and on '
      'the command line `--set section.key=value` does the same (the value is read as YAML). generate() writes the '
      'spec that it built to `<out>/<name>_spec.yaml`.')
    w('')
    w('Named processes (`process=`):')
    w('')
    w('| name | depth h | gap G | post along flow g_p | post across flow l_g | floor | source |')
    w('|---|---|---|---|---|---|---|')
    for n, p in api.PROCESSES.items():
        w('| `%s` | %g | %g | %g | %g | %g | %s |' % (n, p['depth_um'], p['gap_um'], p['post_length_um'], p['post_width_um'],
                                                     p['min_feature_um'], PROCESS_SOURCES[n]))
    w('')
    w('## 2. Part styles at a glance')
    w('')
    w('The tables give each enumerated option, the value to write, what it draws, where it comes from, an example spec '
      'in this repository, and its restrictions. The record field in section 3 verifies the choice on the drawn mask.')
    for (sec, key), choices in SP.CHOICES.items():
        arg = SECTIONS[sec][0]
        w('')
        w('### `%s.%s`  (generate: `%s`)' % (sec, key, arg))
        w('')
        w('| value | draws / does | origin | example | restrictions |')
        w('|---|---|---|---|---|')
        for v in choices:
            d = VALUES[(sec, key, v)]
            w('| `%s` | %s | %s | %s | %s |' % (v, d[0], d[1] or '', ('`%s`' % d[2]) if d[2] else '', d[3] or ''))
    w('')
    w('## 3. Every key')
    w('')
    w('The types and defaults are the code\'s. "Verified by" names the record check (`assessment.checks`) or field '
      '(`validation.*`) that tests the choice on the drawn mask. The symbols follow the frozen glossary. Note that '
      '`post_length_um` is the glossary\'s g_p, "post width", along the flow, and `post_width_um` is its l_g, "link '
      'length", across the flow. Captions, documents, and the dissertation use the glossary\'s words (g_p post width '
      'along the flow, ℓ_g link length across it, Λ row pitch), whereas the YAML keys are historical.')
    for sec, sub in SP.SCHEMA.items():
        arg, what = SECTIONS[sec]
        w('')
        w('### `%s`: %s  (generate: `%s`)' % (sec, what, arg))
        if not isinstance(sub, dict):
            continue
        w('')
        w('| key | type | allowed | default | meaning | symbol | verified by |')
        w('|---|---|---|---|---|---|---|')
        for k, t in sub.items():
            m, sym, ver = KEYS[(sec, k)]
            ch = SP.CHOICES.get((sec, k))
            w('| `%s` | %s | %s | %s | %s | %s | %s |' % (k, _fmt(t), ', '.join('`%s`' % c for c in ch) if ch else '',
                                                       default_of(sec, k), m, sym, ver))
    w('')
    w('## 4. Combination rules (refused with a SpecError, exit 6)')
    w('')
    for r in RULES:
        w('- ' + r)
    w('')
    w('## 5. Example file names')
    w('')
    w('`examples/matrix/mNN_<layout>_<taper>_<head>_<inlet>_<outlet>[_suffix].yaml`:')
    w('')
    w('| field | codes | key |')
    w('|---|---|---|')
    for f, c, k in MATRIX_CODES:
        w('| %s | %s | `%s` |' % (f, c, k))
    w('')
    w('## 6. Exit codes')
    w('')
    w('| code | meaning |')
    w('|---|---|')
    for c, m in EXIT_CODES:
        w('| %s | %s |' % (c, m))
    return '\n'.join(L) + '\n'


RULES = [
    'targets: `cutoff_um` > 0; `filtrate_fraction` in (0, 1), or `filtrate_to_retentate_ratio` > 0 (not both, unless equal); '
    '`f_gap` is never a target (a `targets.f_gap` key is refused)',
    'the lamina delta* = D/kappa_g must be less than half of `width_in_um` and of `width_min_um` (App. A, implementers step 3)',
    'a design whose closure has no solution (f_gap far above the calibrated range) is a spec error rather than a crash',
    '`layout.legs`: 1 for linear; >= 2 for paperclip; odd and >= 3 for serpentine; or `auto`',
    '`outlet.style: ring` only with linear and paperclip layouts',
    '`inlet.style: dinh_filter` only with `meander.style: arms`, `arms: 3`',
    '`meander.seal: allow` is not accepted in strict mode',
    '`start.style: floor` only without a meander head; `start.side_width_um` only with `floor`',
    '`centre_channel.width_min_um` only with `constant_shear` (and required there)',
    '`meander.rectilinear_rows` >= 1',
    '`process.post_width_um` (across the flow) >= `post_length_um` (along it)',
    '`checks.cutoff_tolerance` in (0, 0.05]; `checks.F_tolerance` in (0, 0.01]',
    '`physics.bend_model: none` not in strict mode',
    'a `reproduction` block, `meander.transition: paper` and the other reproduction-only options need `gates.mode: reproduce`',
    'a sealed start without a head must not draw a first side channel below the floor (use a head, `start: floor`, or a finer process)',
    'unknown keys and unknown values are refused',
]
EXIT_CODES = [
    (0, 'PASS (new design) or REPRODUCED (preset)'),
    (2, 'FAIL: a hard check failed on the drawn mask (in strict mode this includes the aspect ratio of the narrowest '
        'CIF channel, whereas the chip area is reported and is never a failure). Since rc7, the photomask is kept for '
        'diagnosis as `<name>_FAILED.dxf` and never under the deliverable name'),
    (3, 'outside the validity envelope: in evaluate or reproduce mode, where the design is drawn and measured '
        '(build), or planned and outlined in memory with nothing written (`opencif check`, `generate(check_only=True)` '
        '/ `generate --check`); and in strict mode when a gate fails only on the drawn mask after the plan passed it '
        '(docs/LIMITATIONS.md item 31)'),
    (4, 'INCOMPLETE (for example, the DRC was not run)'),
    (5, 'REFUSED, with no mask written. Strict mode refuses before drawing outside the envelope or on a failed '
        'pre-flight gate (a meander head\'s GS/G outside gamma_GS [1.4, 1.5], or a channel below w_min). Any mode '
        'refuses a head that cannot seal gap 1 under `seal: require`, an outlet that cannot be planned at any leg '
        'count that the layout tries (found before drawing, or while fitting or drawing), or, since rc7, parts of the '
        'drawn outline that collide (a part closer than the floor to one that it is not joined to, outlet ports not '
        'separated, or an outlet port within 2 mm of the inlet port). The collisions are found on the outline in '
        'memory before the photomask is written, and they are reported with the parts, the remedies, and a figure. The '
        'record gives the numbers and the cutoff band, which on the plateau is the kappa_g interval x delta* and below '
        'it App. B C7\'s plateau-anchored band rather than the self-consistent cutoff, which can lie lower (C7 CQ3)'),
    (6, 'SPEC ERROR (command line) for an invalid or unplannable spec, for which generate() raises '
        'opencif.spec.SpecError'),
    (1, 'an unexpected exception (no record written)'),
]

if __name__ == '__main__':
    open(OUT, 'w').write(render())
    print('wrote %s' % OUT)
