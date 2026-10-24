"""opencif.spec: the user-facing device specification (a YAML file or a dict) and its mapping onto the engine.

The spec has two decoupled targets, the size cutoff D (cutoff_um) and the filtrate fraction F (filtrate_fraction).
D sets the lamina, delta* = D / kappa_g (the cubic kernel pair, kappa_g = 4.9138), so that every gap extracts
f(i) = g(delta*/w_c(i)) and its cutoff is D whatever the width of the center channel. F sets where the march stops,
that is, the number of gaps N. Nothing else is computed by hand. Every option is validated, so an unknown key or value
raises an error, and the values that the build takes from defaults are listed in the record.
"""
import copy
import json
import math

from opencif.physics import schedule as S  # noqa: E402

KAPPA_G = S.KAPPA_G

SCHEMA = {
    'name': str,
    'description': str,
    'targets': {'cutoff_um': float, 'filtrate_fraction': float, 'filtrate_to_retentate_ratio': float},
    'process': {'depth_um': float, 'gap_um': float, 'post_length_um': float, 'post_width_um': float,
                'wall_along_flow_um': float, 'wall_cross_flow_um': float,
                'min_feature_um': float, 'min_radius_um': float, 'fillet_um': float, 'gap_entrance_c': float,
                'aspect_ratio_max': float},
    'chip': {'width_mm': float, 'length_mm': float, 'margin_um': float},
    'centre_channel': {'taper': str, 'width_in_um': float, 'width_min_um': float, 'shear_formula': str},
    'layout': {'style': str, 'legs': (int, str), 'leg_spacing_um': float, 'order': str},
    'meander': {'enabled': bool, 'style': str, 'arms': int, 'channel_um': float, 'wall_um': float,
                'narrow_um': float, 'head_width_um': (float, str), 'seal': str, 'transition': str,
                'rectilinear_rows': int, 'seam_side_width_um': float, 'seam_at': str,
                'turns': str, 'corner_radius_um': float, 'inner_corner_radius_um': float,
                'junction_radius_um': float},
    'inlet': {'style': str, 'port_radius_um': float, 'stem_half_width_um': float, 'port_gap_um': float},
    'outlet': {'style': str, 'port_radius_um': float, 'centre_port_radius_um': float, 'stretch_um': float,
               'flare_deg': float, 'flare_start_um': float,
               'tubing': {'centre_id_mm': float, 'side_id_mm': float, 'length_m': float},
               'ring': {'radius_um': float, 'centre_offset_um': float, 'filtrate_port_offset_um': float,
                        'turn_radius_um': float}},
    'operation': {'flow_uL_min': float},
    'fluid': {'density_kg_m3': float, 'viscosity_Pa_s': float, 'rheology': str},
    'start': {'style': str, 'side_width_um': float},
    'physics': {'fail_safe_margin': bool, 'trim_F': bool, 'wall_rule': str, 'bend_model': str},
    'gates': {'mode': str},
    'checks': {'drc': bool, 'cutoff_tolerance': float, 'F_tolerance': float, 'passage_report': bool},
    'array': {'count': int, 'pitch_um': float, 'inlet': str, 'bus_width_um': float},
    'reproduction': dict,          # engine-level overrides for reproducing a published drawing (flagged in records)
    'paper': dict,                 # free-form source notes and published comparison values (not used by the engine)
}

CHOICES = {
    ('centre_channel', 'taper'): ('constant', 'constant_shear'),
    ('centre_channel', 'shear_formula'): ('exact', 'parallel_plates'),
    ('layout', 'style'): ('linear', 'serpentine', 'paperclip'),
    ('layout', 'order'): ('nested', 'path'),
    ('meander', 'style'): ('pack', 'arms', 'hairpins', 'narrow'),
    ('meander', 'seal'): ('require', 'auto', 'allow'),
    ('meander', 'transition'): ('exact', 'paper'),
    ('meander', 'seam_at'): ('segment', 'post'),
    ('meander', 'turns'): ('arc', 'square'),
    ('inlet', 'style'): ('plain', 'filter', 'dinh_filter'),
    ('outlet', 'style'): ('separate', 'teardrop', 'ring'),
    ('physics', 'wall_rule'): ('exact', 'legacy'),
    ('physics', 'bend_model'): ('stokes', 'none'),
    ('gates', 'mode'): ('strict', 'evaluate', 'reproduce'),
    ('start', 'style'): ('sealed', 'floor'),
    ('array', 'inlet'): ('shared', 'separate'),
    ('fluid', 'rheology'): ('newtonian',),
}

# The glossary's words for the two post dimensions are accepted as aliases of the YAML keys (ADJ-C4). The glossary's
# g_p (post width) is the along-flow extent, post_length_um, and its l_g (link length, the wall thickness that the fluid
# crosses) is the cross-flow extent, post_width_um. load() maps an alias onto its key, and when both are given they must
# agree.
ALIASES = {('process', 'wall_along_flow_um'): 'post_length_um',
           ('process', 'wall_cross_flow_um'): 'post_width_um'}

# REPRO_FIELDS lists the engine fields that a 'reproduction' block may override, and each one is recorded as a
# deviation from the generic part. REPRO_ONLY lists the fields that deliberately depart from the constant-cutoff
# design, and the record lists them as reproduction-only options.
REPRO_FIELDS = {'f_head', 'NM_force', 'R0_arms', 'loop_pull_in', 'loop_pull_in_leg', 'loop_pull_mode', 'loop_ramp',
                'pc_last_short', 'leg_Y', 'inlet_gap', 'pc_inlet_offset', 'pc_out_ext', 'ring_R', 'ring_dc', 'ring_dfp',
                'ring_rq', 'margin', 'n_arms_fixed', 'cap_frac', 'ws_in', 'L_out', 'stem_half', 'wall_M', 'GS',
                'transition', 'ws_rect', 'n_rect', 'seam_at', 'port_r', 'port_r_c', 'port_r_in', 'r_fillet',
                'filter', 'head', 'w_arm', 'w_narrow', 'R_ext_c', 'R_ext_s', 'W', 'depicted_head'}
REPRO_ONLY = {'f_head', 'loop_pull_in', 'pc_last_short', 'leg_Y', 'pc_inlet_offset', 'pc_out_ext', 'NM_force',
              'depicted_head'}


class SpecError(ValueError):
    """The error raised for an invalid specification; it is a ValueError."""
    pass


# CORNER_ROUT_SOLVED is the largest outer-corner rounding r_out / w that the corner table solves (schedule) (DINH-1,
# code review R14).
CORNER_ROUT_SOLVED = S.CORNER_ROUT_SOLVED


def _check(d, schema, path=''):
    """Check the mapping d against the schema at path, raising SpecError for an unknown key or a wrong type and, at the
    top level, for a value outside its CHOICES."""
    if not isinstance(d, dict):
        raise SpecError('%s must be a mapping' % (path or 'spec'))
    for k, v in d.items():
        p = (path + '.' + k) if path else k
        if k not in schema:
            raise SpecError('unknown key %r (allowed: %s)' % (p, ', '.join(sorted(schema))))
        t = schema[k]
        if isinstance(t, dict):
            _check(v, t, p)
        elif t is dict:
            if not isinstance(v, dict):
                raise SpecError('%s must be a mapping' % p)
        else:
            ts = t if isinstance(t, tuple) else (t,)
            if float in ts and isinstance(v, int) and not isinstance(v, bool):
                continue
            if not isinstance(v, ts):
                raise SpecError('%s: expected %s, got %r' % (p, '/'.join(x.__name__ for x in ts), v))
    for (sec, key), allowed in CHOICES.items():
        if path == '' and sec in d and key in d[sec] and d[sec][key] not in allowed:
            raise SpecError('%s.%s = %r is not one of %s' % (sec, key, d[sec][key], allowed))


def load(path_or_dict):
    """Load a spec from a YAML or JSON file (or copy a dict), validate it, and return it as a dict. The checks are the
    schema, the aliases, the fluid, the targets, the required keys, the positive values, and the conflicts with a
    reproduction block, and each failure raises SpecError."""
    if isinstance(path_or_dict, dict):
        d = copy.deepcopy(path_or_dict)
    else:
        txt = open(path_or_dict).read()
        if str(path_or_dict).endswith('.json'):
            d = json.loads(txt)
        else:
            import yaml
            d = yaml.safe_load(txt)
    _check(d, SCHEMA)
    for (sec, alias), key in ALIASES.items():                 # the glossary's aliases (ADJ-C4)
        s_ = d.get(sec, {})
        if alias in s_:
            v = s_.pop(alias)
            if key in s_ and abs(float(s_[key]) - float(v)) > 1e-12:
                raise SpecError('%s.%s (%r) and %s.%s (%r) name the same dimension and disagree: give one'
                                % (sec, alias, v, sec, key, s_[key]))
            s_[key] = v
    fl = d.get('fluid', {})                                   # a positive density and viscosity (ADJ-C5)
    for k_ in ('density_kg_m3', 'viscosity_Pa_s'):
        if k_ in fl and not float(fl[k_]) > 0:
            raise SpecError('fluid.%s must be > 0' % k_)
    t = d.get('targets', {})
    if 'filtrate_to_retentate_ratio' in t:
        # The filtrate-to-retentate flow ratio r = F/(1-F), which is not the glossary's FR (Inglis and Herman's
        # filtration ratio), gives F = r / (1 + r).
        FR = float(t['filtrate_to_retentate_ratio'])
        if not FR > 0:
            raise SpecError('targets.filtrate_to_retentate_ratio must be > 0')
        if 'filtrate_fraction' in t and abs(float(t['filtrate_fraction']) - FR / (1.0 + FR)) > 1e-12:
            raise SpecError('targets.filtrate_fraction and targets.filtrate_to_retentate_ratio disagree: give one')
        t['filtrate_fraction'] = FR / (1.0 + FR)
    for need in (('targets', 'cutoff_um'), ('targets', 'filtrate_fraction'), ('process', 'depth_um'),
                 ('process', 'gap_um'), ('process', 'min_feature_um'), ('centre_channel', 'width_in_um')):
        if need[0] not in d or need[1] not in d[need[0]]:
            raise SpecError('missing required %s.%s' % need)
    D, F = d['targets']['cutoff_um'], d['targets']['filtrate_fraction']
    if not (D > 0) or not math.isfinite(float(D)):          # an infinite D once reached the planner
        raise SpecError('targets.cutoff_um must be a finite number > 0')
    if not (0 < F < 1):
        raise SpecError('targets.filtrate_fraction must lie in (0, 1)')
    positive_values(d)                                        # finite positive values (third-party reviews of rc7)
    conflicts_with_reproduction(d)
    return d


# Every length, depth, flow, radius, tolerance, and count that the design divides by or draws must be a finite
# number > 0 (third-party reviews of rc7, 2026-10-20, review 1 section E; review8's stress test and patch review).
# RC7 accepted a negative or zero operation.flow_uL_min (the Reynolds and Stokes gates then passed on non-positive
# numbers), a zero gap or tubing bore (ZeroDivisionError), a negative depth (refused only inside the planner), a
# negative or zero chip, and an array count of 0 or -2. The values that may be 0 and are refused only when negative
# are process.fillet_um (0 draws square corners, as to_engine's own fillet message advises), process.gap_entrance_c
# (0 gives R_gap without its entrance term, r(G) l_g), the square-turn corner radii and the junction radius (0 gives a
# sharp corner), chip.margin_um, inlet.port_gap_um, outlet.flare_start_um, and the tubing's length (0 means no
# tubing). The ring's two offsets are positions and are not checked here. A count of zero (array.count, meander.arms,
# layout.legs, rectilinear rows) is refused.
POSITIVE = (('process', 'depth_um'), ('process', 'gap_um'), ('process', 'post_length_um'), ('process', 'post_width_um'),
            ('process', 'min_feature_um'), ('process', 'min_radius_um'), ('process', 'aspect_ratio_max'),
            ('chip', 'width_mm'), ('chip', 'length_mm'), ('centre_channel', 'width_in_um'),
            ('centre_channel', 'width_min_um'), ('layout', 'legs'), ('layout', 'leg_spacing_um'),
            ('meander', 'arms'), ('meander', 'channel_um'), ('meander', 'wall_um'), ('meander', 'narrow_um'),
            ('meander', 'head_width_um'), ('meander', 'rectilinear_rows'), ('meander', 'seam_side_width_um'),
            ('inlet', 'port_radius_um'), ('inlet', 'stem_half_width_um'), ('outlet', 'port_radius_um'),
            ('outlet', 'centre_port_radius_um'), ('outlet', 'stretch_um'), ('outlet', 'flare_deg'),
            ('outlet', 'tubing', 'centre_id_mm'), ('outlet', 'tubing', 'side_id_mm'), ('outlet', 'ring', 'radius_um'),
            ('outlet', 'ring', 'turn_radius_um'), ('operation', 'flow_uL_min'), ('start', 'side_width_um'),
            ('checks', 'cutoff_tolerance'), ('checks', 'F_tolerance'), ('array', 'count'), ('array', 'pitch_um'),
            ('array', 'bus_width_um'))
NON_NEGATIVE = (('process', 'fillet_um'), ('process', 'gap_entrance_c'), ('chip', 'margin_um'),
                ('meander', 'corner_radius_um'), ('meander', 'inner_corner_radius_um'),
                ('meander', 'junction_radius_um'), ('inlet', 'port_gap_um'), ('outlet', 'flare_start_um'),
                ('outlet', 'tubing', 'length_m'))


def positive_values(d):
    """Refuse, in one SpecError that lists them all, every value in POSITIVE that is not a finite number > 0 and every
    value in NON_NEGATIVE that is not a finite number >= 0."""
    bad = []
    for keys, strict in ((POSITIVE, True), (NON_NEGATIVE, False)):
        for path in keys:
            v = d
            for k_ in path:
                v = v.get(k_) if isinstance(v, dict) else None
            if v is None or isinstance(v, str):
                continue                                      # absent keys, 'auto' values (the schema checks types)
            if isinstance(v, bool):                           # a flag is not a length or a count (review 3, E9)
                bad.append('%s = %r (must be a number, not true/false)' % ('.'.join(path), v))
                continue
            try:
                fv = float(v)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(fv) or (fv <= 0.0 if strict else fv < 0.0):
                bad.append('%s = %r (must be a finite number %s 0)' % ('.'.join(path), v, '>' if strict else '>='))
    if bad:
        raise SpecError('%s: a length, depth, flow, radius, tolerance or count the design divides by or draws'
                        % '; '.join(bad))


# _ENGINE_TERMS gives the engine value that a spec option's own value stands for, where a reproduction field sets the
# same engine value in other terms (to_engine's mappings; the clearance guarantee, input conflicts).
_ENGINE_TERMS = {('inlet', 'style'): {'plain': False, 'filter': True, 'dinh_filter': 'dinh'},
                 ('meander', 'style'): {'pack': 'pack', 'arms': 'dinh3', 'hairpins': 'dinh3', 'narrow': 'pack'}}


def _same_value(sec, key, given, repro):
    """Return True when a spec value and a reproduction value set the same engine value, with numbers compared to within
    1e-9 of the larger of 1 and their magnitudes."""
    terms = _ENGINE_TERMS.get((sec, key))
    if terms is not None:
        given = terms.get(given, given)
    if isinstance(given, bool) or isinstance(repro, bool) or isinstance(given, str) or isinstance(repro, str):
        return given == repro
    try:
        a_, b_ = float(given), float(repro)
    except (TypeError, ValueError):
        return given == repro
    return abs(a_ - b_) <= 1e-9 * max(1.0, abs(a_), abs(b_))


def conflicts_with_reproduction(d):
    """Refuse a spec that gives a value which its reproduction block replaces with a different one, because the drawing
    would not be what the spec says (the clearance guarantee). Examples are a head width given beside a reproduction W,
    a fillet beside reproduction.r_fillet, and rectilinear rows beside a depicted head's own. Before this check, RC7
    only listed such a value among defaults_used (code review R5b). The remedy is to give only one of the two."""
    bad = []
    for (sec, key) in SPEC_ENGINE:
        if (sec, key) in APPLIES_IF and not APPLIES_IF[(sec, key)](d):
            continue
        if key in d.get(sec, {}):
            o_, v_ = _origin_value(d, sec, key)
            if o_ is not None and not _same_value(sec, key, d[sec][key], v_):
                bad.append('%s.%s = %r is given, but %s replaces it' % (sec, key, d[sec][key], o_))
    repro = d.get('reproduction', {}) or {}
    if 'leg_spacing_um' in d.get('layout', {}) and 'leg_Y' in repro:   # the measured axes place every leg
        bad.append('layout.leg_spacing_um = %r is given, but reproduction.leg_Y = %r places the legs instead'
                   % (d['layout']['leg_spacing_um'], repro['leg_Y']))
    # These are the other spec keys that to_engine maps onto a reproduction field (code review 4, F1).
    pairs = [(('inlet',), 'port_gap_um', 'inlet_gap'), (('outlet',), 'stretch_um', 'L_out'),
             (('outlet', 'ring'), 'radius_um', 'ring_R'), (('outlet', 'ring'), 'centre_offset_um', 'ring_dc'),
             (('outlet', 'ring'), 'filtrate_port_offset_um', 'ring_dfp'), (('outlet', 'ring'), 'turn_radius_um', 'ring_rq')]
    if _m_on(d) and g(d, 'meander', 'style', 'pack') == 'arms':
        pairs.append((('meander',), 'arms', 'n_arms_fixed'))
    for path, key, field in pairs:
        sec = d
        for p_ in path:
            sec = (sec or {}).get(p_) or {}
        if key in sec and field in repro and not _same_value('', key, sec[key], repro[field]):
            bad.append('%s.%s = %r is given, but reproduction.%s = %r replaces it' % ('.'.join(path), key, sec[key], field,
                                                                                  repro[field]))
    # The next two rules cover the two pairs that rc7 LIMITATIONS 29 missed, generalized after review8's patch review
    # (third-party reviews of rc7). The first is an explicit meander.enabled beside a reproduction.head that draws the
    # opposite, because to_engine applies the reproduction block last, so the head drawn was the reproduction's. The
    # second is an explicit width that sets the head's arms beside a different reproduction.w_arm, which replaces the
    # arms' width. That width is meander.narrow_um for the narrow style and meander.channel_um (the head's units, GS)
    # for the pack and arms styles, which is to_engine's w_arm default. A head with arms of their own width stays
    # expressible in the reproduction block alone (reproduction.w_narrow or GS, and w_arm).
    m_ = d.get('meander') if isinstance(d.get('meander'), dict) else {}
    if 'enabled' in m_ and 'head' in repro and not any(b_.startswith('meander.style ') for b_ in bad):
        # The style rule above already names a head that the reproduction replaces, so each conflict gets one message.
        head_on = repro['head'] not in (None, False, '', 'none')
        if bool(m_['enabled']) != head_on:
            bad.append('meander.enabled = %r is given, but reproduction.head = %r %s' % (
                m_['enabled'], repro['head'], 'draws a head' if head_on else 'draws none'))
    st_ = g(d, 'meander', 'style', 'pack')
    k_arm = {'narrow': 'narrow_um', 'pack': 'channel_um', 'arms': 'channel_um'}.get(st_)   # hairpins have no arms
    if (_m_on(d) and k_arm and k_arm in m_ and 'w_arm' in repro
            and not _same_value('', 'w_arm', m_[k_arm], repro['w_arm'])):
        bad.append('meander.%s = %r is given (the %s head draws its arms at it), but reproduction.w_arm = %r replaces '
                   'the arms\' width (give the head\'s widths in the reproduction block, reproduction.%s and w_arm, '
                   'or make them agree)' % (k_arm, m_[k_arm], st_, repro['w_arm'],
                                             'w_narrow' if st_ == 'narrow' else 'GS'))
    tub = g(d, 'outlet', 'tubing', None)
    if tub and ('R_ext_c' in repro or 'R_ext_s' in repro):        # both set the off-chip outlet resistances
        bad.append('outlet.tubing is given, but reproduction.%s sets the off-chip outlet resistance instead'
                   % ' and reproduction.'.join(k for k in ('R_ext_c', 'R_ext_s') if k in repro))
    if bad:
        raise SpecError('%s: the reproduction value would be drawn, not the spec\'s; give one of them' % '; '.join(bad))


def g(d, sec, key, default=None):
    """Return d[sec][key], or default when the section or the key is absent."""
    return d.get(sec, {}).get(key, default)


DEFAULTS = {('process', 'post_length_um'): 25.0, ('process', 'post_width_um'): 50.0,
            ('process', 'gap_entrance_c'): 1.0, ('process', 'aspect_ratio_max'): 10.0, ('chip', 'margin_um'): 500.0, ('centre_channel', 'taper'): 'constant',
            ('centre_channel', 'shear_formula'): 'exact',
            ('meander', 'enabled'): False, ('meander', 'style'): 'pack', ('meander', 'seal'): 'require',
            ('meander', 'transition'): 'exact', ('meander', 'rectilinear_rows'): 2, ('meander', 'seam_at'): 'segment',
            ('meander', 'turns'): 'arc', ('layout', 'order'): 'nested',
            ('inlet', 'style'): 'plain', ('outlet', 'style'): 'separate', ('outlet', 'port_radius_um'): 750.0,
            ('operation', 'flow_uL_min'): 30.0, ('physics', 'fail_safe_margin'): False, ('physics', 'trim_F'): False,
            ('start', 'style'): 'sealed',
            ('checks', 'passage_report'): True,
            ('fluid', 'density_kg_m3'): 1000.0, ('fluid', 'viscosity_Pa_s'): 1.0e-3, ('fluid', 'rheology'): 'newtonian',
            ('physics', 'wall_rule'): 'exact', ('physics', 'bend_model'): 'stokes', ('gates', 'mode'): 'strict', ('layout', 'style'): 'linear'}
DERIVED = {('process', 'min_radius_um'): 'min(4.5, 0.25 x min_feature)',
           ('process', 'fillet_um'): 'min(5, 0.4 x the narrowest drawn feature (min_feature, gap and post length; with a '
                                     'head also the meander wall, channel, arm and narrow widths), 0.3 x the narrowest '
                                     'head unit (channel, arm, narrow width))', ('meander', 'channel_um'): '1.45 x gap',
           ('meander', 'wall_um'): 'max(min_feature, channel/2)', ('meander', 'seam_side_width_um'): 'gap',
           ('meander', 'narrow_um'): 'min_feature', ('meander', 'head_width_um'): 'the narrowest that seals gap 1, searched from 0.3 to 1 x (chip width - 2 margin) (or 3200)',
           ('outlet', 'centre_port_radius_um'): 'port radius', ('inlet', 'port_radius_um'): 'outlet port radius',
           ('inlet', 'stem_half_width_um'): 'w_c/2 + max(gap, min_feature)',
           # the head's smoothing rounds every corner at the fillet at least (code review R1)
           ('meander', 'corner_radius_um'): 'the head fillet (a square turn\'s outer corner is drawn and modelled at '
                                            'max(radius, fillet): the head smoothing rounds it at the fillet)',
           ('meander', 'inner_corner_radius_um'): 'the head fillet (a square turn\'s inner corner is drawn and modelled '
                                                  'at max(radius, fillet): the head smoothing rounds it at the fillet)'}


# A default or derived value is listed among the defaults used only where it applies, so that a record lists what the
# build actually took from a default. This also holds for the meander's values without a head and for the inlet stem
# without a filter or an arm head, which RC6 listed although nothing drew them (code review R5).
def _m_on(d):
    """Return True when the spec enables a meander head."""
    return bool(g(d, 'meander', 'enabled', False))


def _square(d):
    """Return True when the spec's meander head draws square turns."""
    return _m_on(d) and g(d, 'meander', 'turns', 'arc') == 'square'


APPLIES_IF = {('centre_channel', 'shear_formula'): lambda d: g(d, 'centre_channel', 'taper', 'constant') == 'constant_shear',
              ('layout', 'order'): lambda d: g(d, 'layout', 'style', 'linear') == 'paperclip',
              ('meander', 'narrow_um'): lambda d: _m_on(d) and g(d, 'meander', 'style', 'pack') == 'narrow',
              ('meander', 'corner_radius_um'): _square, ('meander', 'inner_corner_radius_um'): _square,
              # the stem beside which the arms rise (schedule Spec.R0), or the filter's half-width (layout)
              ('inlet', 'stem_half_width_um'): lambda d: (g(d, 'inlet', 'style', 'plain') != 'plain' or
                                                          (_m_on(d) and g(d, 'meander', 'style', 'pack') in ('arms', 'hairpins')))}
for _k in ('style', 'seal', 'transition', 'rectilinear_rows', 'seam_at', 'turns', 'channel_um', 'wall_um',
           'seam_side_width_um', 'head_width_um', 'junction_radius_um'):
    APPLIES_IF.setdefault(('meander', _k), _m_on)


# SPEC_ENGINE gives the engine field that each defaulted or derived spec value sets, where a reproduction block can set
# that field instead (spec.REPRO_FIELDS), so that defaults_used names the reproduction value as the one drawn (ADJ-C8).
SPEC_ENGINE = {('process', 'fillet_um'): 'r_fillet', ('meander', 'channel_um'): 'GS', ('meander', 'wall_um'): 'wall_M',
               ('meander', 'seam_side_width_um'): 'ws_in', ('meander', 'narrow_um'): 'w_narrow',
               ('meander', 'head_width_um'): 'W', ('outlet', 'centre_port_radius_um'): 'port_r_c',
               ('inlet', 'port_radius_um'): 'port_r_in', ('inlet', 'stem_half_width_um'): 'stem_half',
               ('meander', 'transition'): 'transition', ('meander', 'rectilinear_rows'): 'n_rect',
               ('meander', 'seam_at'): 'seam_at', ('outlet', 'port_radius_um'): 'port_r', ('inlet', 'style'): 'filter',
               ('meander', 'style'): 'head'}
# DERIVED_ATTR gives the engine attribute that a derived value is read from, for its effective value (ADJ-C8).
DERIVED_ATTR = {('process', 'min_radius_um'): 'rmin_ok', ('process', 'fillet_um'): 'r_fillet',
                ('meander', 'corner_radius_um'): 'r_corner_out', ('meander', 'inner_corner_radius_um'): 'r_corner_in',
                ('meander', 'channel_um'): 'GS', ('meander', 'wall_um'): 'wall_M', ('meander', 'seam_side_width_um'): 'ws_in',
                ('meander', 'narrow_um'): 'w_narrow', ('outlet', 'centre_port_radius_um'): 'port_r_c',
                ('inlet', 'port_radius_um'): 'port_r_in', ('inlet', 'stem_half_width_um'): 'stem_half'}


def _origin_value(d, sec, key):
    """Return (label, value) for the reproduction field that sets this value instead of the spec, its default, or its
    derivation, or (None, None) when no such field does (ADJ-C8)."""
    repro = d.get('reproduction', {}) or {}
    f_ = SPEC_ENGINE.get((sec, key))
    if f_ and f_ in repro:
        return 'reproduction.%s = %r' % (f_, repro[f_]), repro[f_]
    if (sec, key) == ('meander', 'rectilinear_rows') and 'depicted_head' in repro:
        import os
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'depicted.json')) as fh_:
            dj = json.load(fh_)
        n_ = dj.get(repro['depicted_head'], {}).get('n_rect')
        if n_:                                            # to_engine sets n_rect from the depicted head only if > 0
            return 'reproduction.depicted_head = %r (its n_rect %d, data/depicted.json)' % (repro['depicted_head'], n_), n_
    return None, None


def _origin(d, sec, key):
    """Return the label from _origin_value, or None."""
    return _origin_value(d, sec, key)[0]


def defaults_used(d, spec=None):
    """Return every value that the build takes from a default or a derivation rather than from the spec, as listed in
    the record. A value that a reproduction block sets instead is named with that origin, because the reproduction
    value is drawn (ADJ-C8). Given the engine spec, a derived value also carries its effective value, whereas the head
    width searched at build time is in the record's design rather than here."""
    out = []
    for (sec, key), v in DEFAULTS.items():
        if (sec, key) in APPLIES_IF and not APPLIES_IF[(sec, key)](d):
            continue                                      # an option that does not apply to this spec is not "used"
        if key not in d.get(sec, {}):
            o_ = _origin(d, sec, key)
            out.append('%s.%s = %r (default)' % (sec, key, v) if o_ is None else
                       '%s.%s: not the default %r: %s is drawn (rc7)' % (sec, key, v, o_))
    for (sec, key), rule in DERIVED.items():
        if (sec, key) in APPLIES_IF and not APPLIES_IF[(sec, key)](d):
            continue
        if key not in d.get(sec, {}):
            o_ = _origin(d, sec, key)
            if o_ is not None:
                out.append('%s.%s: not derived (%s): %s is drawn (rc7)' % (sec, key, rule, o_))
                continue
            eff = getattr(spec, DERIVED_ATTR[(sec, key)], None) if (spec is not None and (sec, key) in DERIVED_ATTR) else None
            out.append('%s.%s = %s (derived)' % (sec, key, rule) +
                       (' -> %.6g' % float(eff) if isinstance(eff, (int, float)) and not isinstance(eff, bool) else ''))
    # A value that the spec gives and a reproduction block replaces is listed, because the drawn value is the latter. A
    # loaded spec can no longer carry one, because load refuses it (conflicts_with_reproduction), so this lists only
    # specs passed here without load (code review R5b).
    for (sec, key) in SPEC_ENGINE:
        if (sec, key) in APPLIES_IF and not APPLIES_IF[(sec, key)](d):
            continue
        if key in d.get(sec, {}):
            o_, v_ = _origin_value(d, sec, key)
            if o_ is not None and v_ != d[sec][key]:
                out.append('%s.%s = %r given, not drawn: %s is drawn instead (rc7)' % (sec, key, d[sec][key], o_))
    return out


def to_engine(d):
    """Map a validated user spec onto the engine's Spec and return (spec, options)."""
    D = float(d['targets']['cutoff_um'])
    F = float(d['targets']['filtrate_fraction'])
    h = float(d['process']['depth_um'])
    G = float(d['process']['gap_um'])
    g_p = float(g(d, 'process', 'post_length_um', 25.0))
    l_g = float(g(d, 'process', 'post_width_um', 50.0))
    F_min = float(g(d, 'process', 'min_feature_um', min(G, 19.0)))
    wc0 = float(d['centre_channel']['width_in_um'])
    taper = {'constant': 'const', 'constant_shear': 'shear'}[g(d, 'centre_channel', 'taper', 'constant')]
    wc_min = float(g(d, 'centre_channel', 'width_min_um', 0.0))
    if taper == 'shear' and not wc_min:
        raise SpecError('centre_channel.width_min_um is required for the constant_shear taper (its floor)')
    m_on = bool(g(d, 'meander', 'enabled', False))
    style = g(d, 'meander', 'style', 'pack')
    GS = float(g(d, 'meander', 'channel_um', round(1.45 * G, 1)))
    wall_M = float(g(d, 'meander', 'wall_um', max(F_min, 0.5 * GS)))
    chip_w = g(d, 'chip', 'width_mm')
    margin = float(g(d, 'chip', 'margin_um', 500.0))
    hw = g(d, 'meander', 'head_width_um', 'auto')
    W = (float(hw) if hw != 'auto' else (chip_w * 1000.0 - 2 * margin if chip_w else 3200.0))
    head, n_arms, w_narrow = 'none', 3, None
    if m_on:
        head = {'pack': 'pack', 'arms': 'dinh3', 'hairpins': 'dinh3', 'narrow': 'pack'}[style]
        n_arms = {'arms': int(g(d, 'meander', 'arms', 3)), 'hairpins': 0}.get(style, 3)
        if style == 'narrow':
            w_narrow = float(g(d, 'meander', 'narrow_um', F_min))
    inlet = g(d, 'inlet', 'style', 'plain')
    filt = {'plain': False, 'filter': True, 'dinh_filter': 'dinh'}[inlet]
    o_style = g(d, 'outlet', 'style', 'separate')
    outlet = 'separate' if o_style == 'separate' else 'combined'
    outlet_style = 'ring' if o_style == 'ring' else 'teardrop'
    port_r = float(g(d, 'outlet', 'port_radius_um', 750.0))
    port_r_c = float(g(d, 'outlet', 'centre_port_radius_um', port_r))
    port_r_in = float(g(d, 'inlet', 'port_radius_um', port_r))
    stem_half = float(g(d, 'inlet', 'stem_half_width_um', wc0 / 2.0 + max(G, F_min)))
    tub = g(d, 'outlet', 'tubing', None)
    R_ext_c = R_ext_s = 0.0
    if tub:
        Lt = float(tub['length_m'])
        R_ext_c = 128.0 * 1e-3 * Lt / (math.pi * (tub['centre_id_mm'] * 1e-3) ** 4)
        R_ext_s = 128.0 * 1e-3 * Lt / (math.pi * (tub['side_id_mm'] * 1e-3) ** 4)
    kw = dict(name=d.get('name', 'device'), h=h, G=G, g_p=g_p, l_g=l_g, GS=GS, wall_M=wall_M, wc0=wc0,
              delta_star=D / KAPPA_G, F_target=F, ws_in=float(g(d, 'meander', 'seam_side_width_um', G)), W=W,
              F_min=F_min, taper=taper, wc_min=wc_min, head=head, n_arms_fixed=n_arms, w_narrow=w_narrow,
              stem_half=stem_half, filter=filt, outlet=outlet, outlet_style=outlet_style, port_r=port_r,
              port_r_c=port_r_c, port_r_in=port_r_in, r_fillet=float(g(d, 'process', 'fillet_um', 5.0)),
              rmin_ok=float(g(d, 'process', 'min_radius_um', min(4.5, 0.25 * F_min))),
              transition=g(d, 'meander', 'transition', 'exact'), n_rect=int(g(d, 'meander', 'rectilinear_rows', 2)),
              seam_at=g(d, 'meander', 'seam_at', 'segment'), R_ext_c=R_ext_c, R_ext_s=R_ext_s,
              Q_in_uL_min=float(g(d, 'operation', 'flow_uL_min', 30.0)),
              c_ent=float(g(d, 'process', 'gap_entrance_c', 1.0)),
              trim_F=bool(g(d, 'physics', 'trim_F', False)),
              bend_model=str(g(d, 'physics', 'bend_model', 'stokes')),
              shear_formula=str(g(d, 'centre_channel', 'shear_formula', 'exact')),
              turns=str(g(d, 'meander', 'turns', 'arc')),
              r_corner_out=float(g(d, 'meander', 'corner_radius_um', 0.0) or 0.0),
              r_corner_in=float(g(d, 'meander', 'inner_corner_radius_um', 0.0) or 0.0),
              r_junction=float(g(d, 'meander', 'junction_radius_um', 0.0) or 0.0),
              leg_order=str(g(d, 'layout', 'order', 'nested')),
              rho_f=float(g(d, 'fluid', 'density_kg_m3', 1000.0)),
              mu_f=float(g(d, 'fluid', 'viscosity_Pa_s', 1.0e-3)))
    if g(d, 'inlet', 'port_gap_um') is not None:
        kw['inlet_gap'] = float(g(d, 'inlet', 'port_gap_um'))
    if g(d, 'outlet', 'stretch_um') is not None:
        kw['L_out'] = float(g(d, 'outlet', 'stretch_um'))
    if g(d, 'outlet', 'flare_deg') is not None:                   # the flared outlet lead (off by default)
        fd = float(g(d, 'outlet', 'flare_deg'))
        if not 1.0 <= fd <= 45.0:
            raise SpecError('outlet.flare_deg must lie in [1, 45] degrees per side')
        if o_style == 'ring':
            raise SpecError('outlet.flare_deg: the flared lead is drawn for the teardrop and separate outlets only')
        if g(d, 'outlet', 'stretch_um') is not None:              # one source for the stretch's length (rc5 review)
            raise SpecError('outlet.stretch_um with outlet.flare_deg: the flare sets the stretch (outlet.flare_start_um '
                            'straight, then the flare and its end arc to the lead width); give flare_start_um instead')
        if bool(g(d, 'physics', 'fail_safe_margin', False)):
            raise SpecError('outlet.flare_deg with physics.fail_safe_margin is not supported: the fail-safe re-solve '
                            'rebuilds the march, which moves the last side channel\'s width the flare was planned on')
        kw['flare_deg'] = fd
        # The default straight run is Dinh's Fig. 2B read at the corrected scale bar (rc5 read 67).
        kw['flare_start'] = float(g(d, 'outlet', 'flare_start_um', 65.4))
        if kw['flare_start'] <= 0.0:
            raise SpecError('outlet.flare_start_um must be positive')
    elif g(d, 'outlet', 'flare_start_um') is not None:
        raise SpecError('outlet.flare_start_um needs outlet.flare_deg (the straight run before a flare)')
    if outlet_style == 'ring':
        rg = g(d, 'outlet', 'ring', {}) or {}
        # The parametric ring (opencif F6) is scaled on the ports, so that it encloses the center port and its far side
        # passes through the filtrate port, with Xia 2016 Fig. 3a's proportions:
        #     dc  = 0.284 R
        #     dfp = dc + R - 0.353 r_port
        R_def = max(2.2 * port_r_c, port_r_c + port_r + 400.0)
        Rr = float(rg.get('radius_um', R_def))
        kw.update(ring_R=Rr, ring_dc=float(rg.get('centre_offset_um', 0.284 * Rr)),
                  ring_dfp=float(rg.get('filtrate_port_offset_um', 0.284 * Rr + Rr - 0.353 * port_r)),
                  ring_rq=float(rg.get('turn_radius_um', 250.0)))
    repro = d.get('reproduction', {}) or {}
    bad = set(repro) - REPRO_FIELDS
    if bad:
        raise SpecError('reproduction: unknown engine fields %s' % sorted(bad))
    # The flare sets the outlet stretch, so a reproduction.L_out beside it is refused (rc5 verification, deferred).
    if 'L_out' in repro and kw.get('flare_deg') is not None:
        raise SpecError('reproduction.L_out with outlet.flare_deg: the flare sets the outlet stretch '
                        '(outlet.flare_start_um straight, then the flare and its end arc); give flare_start_um instead')
    for k, v in repro.items():
        if k == 'depicted_head':
            # The head is drawn as the paper depicts it, with the prescribed gap fractions measured on its figure
            # (data/depicted.json).
            import os
            dj = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'depicted.json')))[v]
            kw['f_head'] = tuple(dj['f_head'])
            if dj.get('n_rect'):
                kw['NM_force'], kw['n_rect'] = int(dj['NM']), int(dj['n_rect'])
            continue
        kw[k] = tuple(v) if isinstance(v, list) else v
    spec = S.Spec(**kw)
    # In strict mode an outlet stretch-end beyond the solved 2-D tables is refused (rc5 review P, deferred).
    spec.exp_domain_refuse = (g(d, 'gates', 'mode', 'strict') == 'strict')
    for wlab, wv in (('centre_channel.width_in_um', spec.wc0), ('centre_channel.width_min_um', spec.wc_min)):
        if wv and spec.delta_star / wv >= 0.5:
            raise SpecError('the lamina delta* = D/kappa_g = %.3f um is at least half of %s = %.1f um: the cubic skim '
                            'kernel is defined only for delta/w_c < 1/2 (App. A Result 2)' % (spec.delta_star, wlab, wv))
    if spec.head != 'none' and int(spec.n_rect) < 1:          # also when set through a reproduction block
        raise SpecError("rectilinear rows: at least 1 under a meander head's transition unit (n_rect %r)" % spec.n_rect)
    spec.tubing_centre_id_mm = float(tub['centre_id_mm']) if tub else None     # for the tubing scenarios
    if spec.w_arm is None or (w_narrow and 'w_arm' not in repro):
        spec.w_arm = w_narrow if w_narrow else spec.GS          # derived default recomputed (audit F7)
    # The head and outlet smoothing (shapes.smooth, a closing and then an opening of radius r_fillet) fills every solid
    # wall and erases every channel no wider than 2 r_fillet (Reviewer A's E5). The earlier default of 5 um did so at a
    # 10 um floor, where the narrow head's 10 um stems vanished through the post wall and its first units were drawn as
    # islands (E5). The default is now min(5, 0.4 x the narrowest drawn feature, 0.3 x the narrowest head unit), which
    # is still 5 um for every shipped spec that uses the default (the others set 1.5 to 5 um, each below half their
    # floor), and an explicit fillet at or above half the narrowest feature is refused before drawing.
    narrow_ = [('process.min_feature_um', spec.F_min), ('process.gap_um', spec.G), ('process.post_length_um', spec.g_p)]
    if spec.head != 'none':
        narrow_ += [('meander.wall_um', spec.wall_M), ('meander.channel_um', spec.GS), ('the arm width', spec.w_arm)]
        if spec.w_narrow:
            narrow_.append(('meander.narrow_um', spec.w_narrow))
    lab_, lim_ = min(narrow_, key=lambda t: t[1])
    # The head T-junction's fillet term (schedule, S_T) is fitted for r_fillet / w up to 0.3, where w is the unit's
    # width, and it is capped there. The default stays inside that range as well (an independent check of rc6).
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
        notes_spec.append('fillet %.2f um is %.2f of the narrowest head unit (%.1f um): the head T-junction\'s fillet term '
                          '(schedule, S_T) is fitted to r_f / w 0.3 and held there beyond it'
                          % (spec.r_fillet, spec.r_fillet / min(units_), min(units_)))
    lay = g(d, 'layout', 'style', 'linear')
    legs = g(d, 'layout', 'legs', 1 if lay == 'linear' else 'auto')
    if lay == 'linear' and legs not in (1, 'auto'):
        raise SpecError('layout.legs must be 1 for a linear layout')
    if lay == 'serpentine' and isinstance(legs, int) and legs % 2 == 0:
        raise SpecError('layout.legs must be odd for a serpentine (the outlet must leave away from the head)')
    if isinstance(legs, str) and legs != 'auto':
        raise SpecError("layout.legs must be an integer or 'auto'")
    if inlet == 'dinh_filter' and not (m_on and style == 'arms' and n_arms == 3):
        raise SpecError("inlet.style: dinh_filter is Dinh 2024's filter in its own absolute coordinates: it fits only "
                        "Dinh's three-arm head (meander.style: arms, arms: 3)")
    if o_style == 'ring' and g(d, 'layout', 'style', 'linear') == 'serpentine':
        raise SpecError('outlet.style: ring is drawn for linear and paperclip layouts only (its clearance from a '
                        "serpentine's far loop is not checked)")
    mode_ = g(d, 'gates', 'mode', 'strict')
    if mode_ == 'strict' and g(d, 'meander', 'seal') == 'allow':
        raise SpecError("meander.seal: allow (a raised gap-1 fraction, its cutoff above D) is not allowed in strict mode; "
                        "use seal: require, a wider head, the narrow style, or gates.mode: evaluate")
    if mode_ != 'reproduce' and d.get('reproduction'):
        raise SpecError('a reproduction block (engine-level overrides for drawing a published device) needs gates.mode: reproduce')
    tol_ = g(d, 'checks', 'cutoff_tolerance', 0.01)
    if not (0 < tol_ <= 0.05):
        raise SpecError('checks.cutoff_tolerance must lie in (0, 0.05]: a looser cutoff tolerance does not validate D')
    lg_ = g(d, 'layout', 'legs', None)
    if isinstance(lg_, int) and (lg_ < 1 or (g(d, 'layout', 'style') == 'paperclip' and lg_ < 2)
                                 or (g(d, 'layout', 'style') == 'serpentine' and lg_ < 3)):
        raise SpecError('layout.legs: at least 1 (linear), 2 (paperclip) or 3 (serpentine, odd)')
    pl_, pw_ = g(d, 'process', 'post_length_um', 25.0), g(d, 'process', 'post_width_um', 50.0)
    if pw_ < pl_:
        raise SpecError('process.post_width_um (%.1f, across the flow) must be >= post_length_um (%.1f, along the flow): '
                        'posts are stadiums elongated across the flow, and the gap model assumes straight gap walls '
                        'post_width_um long' % (pw_, pl_))
    start = g(d, 'start', 'style', 'sealed')
    if start == 'floor' and m_on:
        raise SpecError("start.style: floor is the no-meander start rule (Gifford 2014); a meander head seals the start itself")
    if m_on and int(g(d, 'meander', 'rectilinear_rows', 2)) < 1:
        raise SpecError("meander.rectilinear_rows must be >= 1: the transition unit's drop runs over the rectilinear "
                        "rows' wall blocks; with none, the straight rows' post band would begin under the unit's quarter "
                        "arc (not drawn or validated)")
    if start != 'floor' and g(d, 'start', 'side_width_um') is not None:
        raise SpecError('start.side_width_um applies only to start.style: floor')
    if g(d, 'centre_channel', 'taper', 'constant') == 'constant' and g(d, 'centre_channel', 'shear_formula') is not None:
        raise SpecError('centre_channel.shear_formula applies only to the constant_shear taper (rc7)')
    if g(d, 'centre_channel', 'shear_formula', 'exact') == 'parallel_plates' and g(d, 'gates', 'mode', 'strict') != 'reproduce':
        raise SpecError('centre_channel.shear_formula: parallel_plates (the parallel-plate shear 6 Q/(a b^2)) departs from '
                        'spec v5 S16-b ("max gamma_dot via C(alpha) [NEVER 6]"): it reads a published device (Xia 2016, '
                        'register XIA-1) and needs gates.mode: reproduce (rc7)')
    if g(d, 'layout', 'order') is not None and g(d, 'layout', 'style', 'linear') != 'paperclip':
        raise SpecError('layout.order applies only to the paperclip layout (rc7)')
    if g(d, 'meander', 'turns', 'arc') == 'square' and not m_on:
        raise SpecError('meander.turns: square applies to a meander head (meander.enabled: true) (rc7)')
    if g(d, 'meander', 'turns', 'arc') != 'square' and (g(d, 'meander', 'corner_radius_um') is not None
                                                       or g(d, 'meander', 'inner_corner_radius_um') is not None):
        raise SpecError('meander.corner_radius_um and inner_corner_radius_um apply only to meander.turns: square (rc7)')
    rj_ = g(d, 'meander', 'junction_radius_um')
    if rj_ is not None:
        if not m_on:
            raise SpecError('meander.junction_radius_um applies to a meander head (rc7)')
        lim_j = 0.5 * min(float(g(d, 'process', 'post_width_um', 50.0)),
                          float(g(d, 'meander', 'wall_um', max(F_min, 0.5 * GS))))
        if not (0.0 <= float(rj_) <= lim_j):
            raise SpecError('meander.junction_radius_um must lie in [0, %.1f] um: half the thinner of the wall the head '
                            'gaps cross (process.post_width_um) and the wall between units (meander.wall_um) (rc7)' % lim_j)
    for k_ in ('corner_radius_um', 'inner_corner_radius_um'):
        v_ = g(d, 'meander', k_)
        if v_ is not None and not (math.isfinite(float(v_)) and float(v_) >= 0.0):
            raise SpecError('meander.%s must be a finite radius >= 0 um (rc7)' % k_)
    if spec.turns == 'square' and units_:
        # The head zone's smoothing (shapes.smooth, a closing and then an opening of radius r_fillet) rounds every
        # corner of the head at r_fillet at least, so a square turn's corners are drawn, measured
        # (shapes.corner_round_area), and modeled (schedule.c_corner) at max(radius, r_fillet) (code review R1).
        ro_, ri_ = spec.r_corner_out, spec.r_corner_in
        spec.r_corner_out, spec.r_corner_in = max(ro_, spec.r_fillet), max(ri_, spec.r_fillet)
        for lab_r, given_, eff_ in (('corner_radius_um', ro_, spec.r_corner_out),
                                    ('inner_corner_radius_um', ri_, spec.r_corner_in)):
            if eff_ > given_ and g(d, 'meander', lab_r) is not None:
                notes_spec.append('meander.%s %.2f um is below the head fillet %.2f um: the head smoothing rounds the '
                                  'corner at the fillet, so it is drawn and modelled at %.2f um (rc7)'
                                  % (lab_r, given_, spec.r_fillet, eff_))
        # The radii are checked against the head's drawn unit widths rather than channel_um's default. The outer
        # rounding is solved to r_out / w 0.307 (provenance/rc7_scripts/corner/corner_table.json), and an inner
        # rounding must fit the wall it rounds (two at a wall's end) and half the unit (code review R14).
        w_sq = min(units_)
        src_f = 'reproduction.r_fillet' if 'r_fillet' in repro else 'process.fillet_um'
        src_o = 'meander.corner_radius_um' if spec.r_corner_out == ro_ else 'the head fillet (%s)' % src_f
        if spec.r_corner_out > CORNER_ROUT_SOLVED * w_sq + 1e-9:
            raise SpecError('%s %.2f um is %.3f of the narrowest head unit (%.2f um): the square corner\'s outer rounding '
                            'is solved to r / w %.3f only (rc7)' % (src_o, spec.r_corner_out, spec.r_corner_out / w_sq,
                                                                    w_sq, CORNER_ROUT_SOLVED))
        lim_i = 0.5 * min(w_sq, spec.wall_M)
        if spec.r_corner_in > lim_i + 1e-9:
            src_i = 'meander.inner_corner_radius_um' if spec.r_corner_in == ri_ else 'the head fillet (%s)' % src_f
            raise SpecError('%s %.2f um exceeds %.2f um: half the narrower of the narrowest head unit (%.2f um) and the '
                            'wall between units (%.2f um; two inner roundings at a wall\'s end) (rc7)'
                            % (src_i, spec.r_corner_in, lim_i, w_sq, spec.wall_M))
    if g(d, 'centre_channel', 'taper', 'constant') == 'constant' and g(d, 'centre_channel', 'width_min_um') is not None:
        raise SpecError('centre_channel.width_min_um applies only to the constant_shear taper (a constant taper keeps '
                        'width_in_um throughout)')
    ar_ = g(d, 'process', 'aspect_ratio_max', 10.0)
    if not (isinstance(ar_, (int, float)) and math.isfinite(ar_) and ar_ > 0):
        raise SpecError('process.aspect_ratio_max must be a finite number > 0 (the deepest channel the process makes, as '
                        'depth over width)')
    F_tol_ = g(d, 'checks', 'F_tolerance', 1e-3)
    if not (0 < F_tol_ <= 0.01):
        raise SpecError('checks.F_tolerance must lie in (0, 0.01]: a looser tolerance does not validate F')
    if g(d, 'physics', 'bend_model', 'stokes') == 'none' and mode_ == 'strict':
        raise SpecError("physics.bend_model: none (the plain 1-D rule at bends) is for comparison only; not in strict mode")
    mode = g(d, 'gates', 'mode', 'strict')
    repro_only = sorted(set(repro) & REPRO_ONLY)
    if spec.transition == 'paper':
        repro_only.append('transition: paper')
    if spec.shear_formula == 'parallel_plates':                   # a departure from spec v5 S16-b (XIA-1)
        repro_only.append('shear_formula: parallel_plates')
    if repro_only and mode != 'reproduce':
        raise SpecError('reproduction-only options %s need gates.mode: reproduce (they deliberately depart from the '
                        'constant-cutoff design)' % repro_only)
    st_ = (d.get('paper', {}) or {}).get('source_trace')          # a pointer to the trace record, if any (ADJ-C6)
    if isinstance(st_, list) and st_ and all(isinstance(x_, str) and x_.strip() for x_ in st_):
        st_ = '; '.join(x_.strip() for x_ in st_)
    if st_ is not None and not (isinstance(st_, str) and st_.strip()):
        # The paper block is free-form (RC6 accepted anything there), so a value that is not a pointer is noted rather
        # than refused (R16).
        notes_spec.append('paper.source_trace %r is not a string or a list of strings: not taken as a trace pointer (rc7)'
                          % (st_,))
        st_ = None
    # The targets that a reproduction fitted to the paper's figures are free text in paper.*. They are copied to
    # assessment.status.reproduction, which otherwise says "D and F from printed values" (codereview3 follow-up).
    tf_ = (d.get('paper', {}) or {}).get('targets_fitted')
    if isinstance(tf_, list) and tf_ and all(isinstance(x_, str) and x_.strip() for x_ in tf_):
        tf_ = '; '.join(x_.strip() for x_ in tf_)
    if tf_ is not None and not (isinstance(tf_, str) and tf_.strip()):
        notes_spec.append('paper.targets_fitted %r is not a string or a list of strings: not taken (rc7)' % (tf_,))
        tf_ = None
    opts = dict(layout=lay, legs=legs, leg_spacing=g(d, 'layout', 'leg_spacing_um'),
                seal=g(d, 'meander', 'seal', 'require'), mode=mode, failsafe=bool(g(d, 'physics', 'fail_safe_margin', False)),
                wall_rule=g(d, 'physics', 'wall_rule', 'exact'), drc=bool(g(d, 'checks', 'drc', True)),
                cut_tol=float(g(d, 'checks', 'cutoff_tolerance', 0.01)), F_tol=float(g(d, 'checks', 'F_tolerance', 1e-3)),
                chip=(g(d, 'chip', 'width_mm'), g(d, 'chip', 'length_mm')), margin=margin,
                head_width_auto=(hw == 'auto' and 'W' not in (d.get('reproduction', {}) or {})), style=(style if m_on else None), repro_only=repro_only,
                inlet=inlet, outlet=o_style, defaults_used=defaults_used(d, spec), start=start,
                start_width=float(g(d, 'start', 'side_width_um', F_min)),
                AR_max=float(ar_), w_min=max(F_min, spec.h / float(ar_)), notes_spec=notes_spec,
                source_trace=st_, targets_fitted=tf_, passage_report=bool(g(d, 'checks', 'passage_report', True)))
    return spec, opts
