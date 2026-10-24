"""rc6: the hashed-code edits other than the rc5 verification's deferred ones (apply_rc5_deferred.py) and the
devspecs.py replica values, in place, each by exact string replacement (it stops if a target is not found exactly
once).  Re-applied on 2026-10-17 after the first rc6 working tree was lost (the edits as recorded in the rc6
reconstruction record), on the verified rc5 release (code 56e1980d8a064b35).

  E5 (Reviewer A): measure2's width guard (a planned feature the drawing does not contain is a named FAIL, not a crash),
     build.py's validation.features_found, pipeline.assess's hard check; the fillet default and its refusal (spec.py)
  the Dinh flare onset default 67 -> 65.4 um; stale comments; the Xia depicted head at the rc6 scale (depicted.py)
  the named processes (api.PROCESSES) at the rc6 values; __version__ 1.0.0rc6

Usage: python3 apply_rc6_code.py OPENCIF_ROOT"""
import os
import sys

R = sys.argv[1]


def edit(rel, old, new, count=1):
    p = os.path.join(R, rel)
    s = open(p).read()
    assert s.count(old) == count, (rel, old[:70], s.count(old))
    open(p, 'w').write(s.replace(old, new))


edit('opencif/__init__.py', "__version__ = '1.0.0rc5'", "__version__ = '1.0.0rc6'")

# ---- api.py: the named processes at the rc6 values (data/r6/rc6_values.json)
edit('opencif/api.py', "    'dinh2024': dict(depth_um=140.0, gap_um=19.0, post_length_um=27.0, post_width_um=54.0, min_feature_um=19.0),",
     "    'dinh2024': dict(depth_um=140.0, gap_um=19.0, post_length_um=26.1, post_width_um=54.0, min_feature_um=19.0),")
edit('opencif/api.py', "    'xia2016': dict(depth_um=140.0, gap_um=16.0, post_length_um=40.25, post_width_um=78.1, min_feature_um=16.0,",
     "    'xia2016': dict(depth_um=140.0, gap_um=16.0, post_length_um=35.8, post_width_um=70.9, min_feature_um=16.0,")
edit('opencif/api.py', "    'gifford2014': dict(depth_um=80.0, gap_um=20.5, post_length_um=25.6, post_width_um=56.8, min_feature_um=19.0),",
     "    'gifford2014': dict(depth_um=80.0, gap_um=20.5, post_length_um=25.3, post_width_um=59.5, min_feature_um=19.0),")

# ---- build.py: the features the measurement did not find
edit('opencif/build.py', """    m = M2.measure_device(fluid, sc_f, info)
    nw = M2.network(sc_f, m)""", """    m = M2.measure_device(fluid, sc_f, info)
    # rc6 (E5): features the measurement probes did not find where the plan puts them (measure2.found)
    V['features_found'] = dict(ok=not m['not_found'], not_found=m['not_found'][:20], n_not_found=len(m['not_found']))
    if m['not_found']:
        f0 = m['not_found'][0]
        say('[%s] DRAWING DEFECT: %d planned feature(s) not found on the drawing (first: row %d, %s, side %s)'
            % (name, len(m['not_found']), f0['row'], f0['kind'], f0['side']))
    nw = M2.network(sc_f, m)""")

# ---- pipeline.py: the hard check; the spec's notes reach the record; outlet_precheck's docstring
edit('opencif/pipeline.py', """    fid = V['drawing_fidelity']
    f_no = fid.get('f_dev_max_where_no_cutoff', fid['f_dev_over_fstar_max_where_no_cutoff'])""",
     """    uf = V.get('features_found')
    if uf is not None:                                   # rc6 (E5): a planned feature the drawing does not contain
        chk('every planned feature found on the drawing', uf['ok'],
            dict(n_not_found=uf['n_not_found'], first=uf['not_found'][:3]) if not uf['ok'] else 'all',
            'every measurement section (head units, gaps, side segments, centre, loops, outlet) finds fluid where '
            'the plan puts the feature')
    fid = V['drawing_fidelity']
    f_no = fid.get('f_dev_max_where_no_cutoff', fid['f_dev_over_fstar_max_where_no_cutoff'])""")
edit('opencif/pipeline.py', """def _prepare(spec, opts, say):
    notes = []""", """def _prepare(spec, opts, say):
    notes = list(opts.get('notes_spec', []))             # rc6: what spec.to_engine found worth saying (none refused)""")
edit('opencif/pipeline.py', """    generate(check_only=True) also report.  A count that fails for another reason (too short for n legs) ends the""",
     """    generate(check_only=True) also report.  rc6: in strict mode an outlet whose stretch-end term lies beyond the solved
    2-D tables is refused here too (layout.OutletOutsideTables); the other modes build it and say so in the verdict.  A
    count that fails for another reason (too short for n legs) ends the""")

# ---- spec.py: the fillet default (E5's cause) and its refusal; the flare onset default
edit('opencif/spec.py', """DEFAULTS = {('process', 'post_length_um'): 25.0, ('process', 'post_width_um'): 50.0, ('process', 'fillet_um'): 5.0,
""", """DEFAULTS = {('process', 'post_length_um'): 25.0, ('process', 'post_width_um'): 50.0,
""")
edit('opencif/spec.py', """DERIVED = {('process', 'min_radius_um'): 'min(4.5, 0.25 x min_feature)', ('meander', 'channel_um'): '1.45 x gap',
""", """DERIVED = {('process', 'min_radius_um'): 'min(4.5, 0.25 x min_feature)',
           ('process', 'fillet_um'): 'min(5, 0.4 x the narrowest drawn feature (min_feature, gap and post length; with a '
                                     'head also the meander wall, channel, arm and narrow widths), 0.3 x the narrowest '
                                     'head unit (channel, arm, narrow width))', ('meander', 'channel_um'): '1.45 x gap',
""")
edit('opencif/spec.py', """        kw['flare_start'] = float(g(d, 'outlet', 'flare_start_um', 67.0))""",
     """        kw['flare_start'] = float(g(d, 'outlet', 'flare_start_um', 65.4))     # rc6: Dinh Fig. 2B at the corrected bar (rc5 67)""")
edit('opencif/spec.py', """    if spec.w_arm is None or (w_narrow and 'w_arm' not in repro):
        spec.w_arm = w_narrow if w_narrow else spec.GS          # derived default recomputed (audit F7)
""", """    if spec.w_arm is None or (w_narrow and 'w_arm' not in repro):
        spec.w_arm = w_narrow if w_narrow else spec.GS          # derived default recomputed (audit F7)
    # rc6 (Reviewer A's E5): the head and outlet smoothing (geometry.smooth: a closing then an opening of radius
    # r_fillet) fills every solid wall and erases every channel no wider than 2 r_fillet.  The default 5 um did so at a
    # 10 um floor (E5: the narrow head's 10 um stems vanished through the post wall and its first units were drawn as
    # islands).  The default is now min(5, 0.4 x the narrowest drawn feature, 0.3 x the narrowest head unit) -- still
    # 5 um for every shipped spec that uses the default (the others set 1.5-5 um, each below half their floor) -- and
    # an explicit fillet at or above half the narrowest feature is refused before drawing
    narrow_ = [('process.min_feature_um', spec.F_min), ('process.gap_um', spec.G), ('process.post_length_um', spec.g_p)]
    if spec.head != 'none':
        narrow_ += [('meander.wall_um', spec.wall_M), ('meander.channel_um', spec.GS), ('the arm width', spec.w_arm)]
        if spec.w_narrow:
            narrow_.append(('meander.narrow_um', spec.w_narrow))
    lab_, lim_ = min(narrow_, key=lambda t: t[1])
    # the head T-junction's fillet term (cifsched, S_T) is fitted for r_fillet / w up to 0.3 (and capped there), w the
    # unit's width: the default stays inside that range too (an independent check of rc6)
    units_ = [w_ for w_ in ((spec.GS, spec.w_arm, spec.w_narrow) if spec.head != 'none' else ()) if w_]
    notes_spec = []
    if g(d, 'process', 'fillet_um') is None and 'r_fillet' not in repro:
        spec.r_fillet = min([5.0, 0.4 * lim_] + [0.3 * w_ for w_ in units_])
    elif spec.r_fillet > 0.0 and spec.r_fillet >= 0.5 * lim_:
        src_ = 'reproduction.r_fillet' if 'r_fillet' in repro else 'process.fillet_um'
        raise SpecError('%s %.2f um is not below half the narrowest drawn feature (%s, %.2f um): the head and outlet '
                        'smoothing (a closing then an opening of that radius) would fill walls and erase channels no '
                        'wider than %.2f um.  Use a fillet < %.2f um (0 draws square corners)'
                        % (src_, spec.r_fillet, lab_, lim_, 2.0 * spec.r_fillet, 0.5 * lim_))
    elif units_ and spec.r_fillet > 0.3 * min(units_):
        notes_spec.append('fillet %.2f um is %.2f of the narrowest head unit (%.1f um): the head T-junction\\'s fillet term '
                          '(cifsched, S_T) is fitted to r_f / w 0.3 and held there beyond it'
                          % (spec.r_fillet, spec.r_fillet / min(units_), min(units_)))
""")
edit('opencif/spec.py', """                AR_max=float(ar_), w_min=max(F_min, spec.h / float(ar_)))""",
     """                AR_max=float(ar_), w_min=max(F_min, spec.h / float(ar_)), notes_spec=notes_spec)""")

# ---- cifsched.py: the flare onset default; a stale comment
edit('opencif/engine/cifsched.py', """    flare_start: float = 67.0        # straight run (um) from the last gap's centre to the flare's onset (Dinh 2024
                                     # Fig. 2B, fourth panel: 67 +- 2 um, provenance/rc5_scripts/measure_dinh_flare.py)""",
     """    flare_start: float = 65.4        # straight run (um) from the last gap's centre to the flare's onset (Dinh 2024
                                     # Fig. 2B, fourth panel: 65.4 +- 2 um at the bar read with its anti-aliased ends,
                                     # provenance/rc6_scripts/measure_rc6.py; rc5 read 67 at the short bar)""")
edit('opencif/engine/cifsched.py', """        # g_p + G/2 long (with Dinh's GS 28, g_p 27, G 19: 50.5 um; measured on Dinh's panel: 49 um)""",
     """        # g_p + G/2 long (with Dinh's GS 28, g_p 26.1, G 19: 49.6 um, rc6; rc5's g_p 27 gave 50.5; measured on Dinh's
        # panel: 49 um at rc5's short bar reading, 48 um at the bar read with its anti-aliased ends)""")

# ---- depicted.py: the Xia head at the rc6 scale and Fig. 3a's corrected lengths
edit('opencif/engine/depicted.py', """micrograph (r4/meas_xia_head2.py: 10 rows at 1.953 um/px; rows 11-12 extrapolated, Fig. 1b(i) shows 12).\"\"\"""",
     """micrograph (r4/meas_xia_head2.py: 10 rows; rc6: at 1.798 um/px, the scale panels b-g's printed widths and loop
curvature give -- data/r6/rc6_values.json; rc5 read them at the panel-i bar, 1.953 um/px).\"\"\"""")
edit('opencif/engine/depicted.py', """# ---- Xia 2016 Fig. 3b micrograph (native px; axis 242; 1.953 um/px) ----------------------------------------------
XIA_SC, XIA_AX = 250.0 / 128.0, 242.0""", """# ---- Xia 2016 Fig. 3b micrograph (native px; axis 242; rc6: 1.798 um/px, data/r6/rc6_values.json) ----------------
_R6 = json.load(open(os.path.join(HERE, '..', 'data', 'r6', 'rc6_values.json')))
XIA_SC, XIA_AX = _R6['xia']['scale_um_per_px'], 242.0""")
edit('opencif/engine/depicted.py', """XIA_HEAD_TOP_MM, XIA_FIRST_PILL_MM = _M3A['head_top_below_inlet_mm'], _M3A['first_pill_below_inlet_mm']   # 1.45, 13.73""",
     """# rc6: Fig. 3a's bar read with its anti-aliased ends (178.32 px; r4 used 179): every Fig. 3a length x 1.003827
XIA_HEAD_TOP_MM = _R6['xia']['layout']['head_top_below_inlet_mm']          # 1.453 (rc5 1.447)
XIA_FIRST_PILL_MM = _R6['xia']['layout']['first_pill_below_inlet_mm']      # 13.786 (rc5 13.734)""")
print('applied to', R)
