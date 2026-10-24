"""opencif 1.0.0: the physics auditor (opencif/checks/auditor.py) in the build, the assessment, the record and the CLI.
Applied after refactor_v1.py and patch_v1.py.  The photomask is not touched: every edit is after the DXF is written.
Blocks of build() that the standalone audit (`opencif audit RUN_DIR`) needs as well are moved, verbatim, into
functions of their own (rounding_check, outlet_and_loops_measured, drawn_parameters, intended_and_declared, run_drc).
usage: python -I patch_audit.py TREE AUDITOR_PY"""
import os
import shutil
import sys

T, AUD = sys.argv[1:3]
n_done = 0


def edit(rel, old, new, count=1):
    global n_done
    p = os.path.join(T, rel)
    s = open(p, encoding='utf-8').read()
    assert s.count(old) == count, (rel, s.count(old), old[:100])
    s = s.replace(old, new)
    open(p, 'w', encoding='utf-8').write(s)
    n_done += 1


shutil.copy(AUD, os.path.join(T, 'opencif', 'checks', 'auditor.py'))

# ------------------------------------------------------------------------------------------------------------- build.py
B = 'opencif/build.py'
edit(B, "from opencif.checks import measure as M2  # noqa: E402\n",
     "from opencif.checks import measure as M2  # noqa: E402\nfrom opencif.checks import auditor as AU  # noqa: E402\n")
edit(B, '''# CIF_v2's embedded physics auditor (App. B "P8"; CIF_v2 _PHYSICS_SRC, PhysicsSpec.Re_max, .Q_in_uLmin and check
# P8_stokes_gates): its inlet channel limit, its default operating point and its entrance-length criterion
# 0.05 Re D_h <= Lambda/2.  Reported in every record at both flows, not gated (docs/LIMITATIONS.md 6);
# tests/test_appendix_trace.py pins the values to the CIF_v2 source text.
P8_RE_MAX = 1.0
P8_Q_UL_MIN = 5.0''', '''# CIF_v2's embedded physics auditor (App. B "P8"; CIF_v2 _PHYSICS_SRC, PhysicsSpec.Re_max, .Q_in_uLmin and check
# P8_stokes_gates): its own inlet channel limit, its default operating point and its entrance-length criterion
# 0.05 Re D_h <= Lambda/2.  The port (opencif.checks.auditor, check P8) evaluates them on the drawn mask; gates() gives
# their design-stage values.  Reported, not gated (docs/LIMITATIONS.md 6); tests/test_appendix_trace.py pins the values
# to the CIF_v2 source text.
P8_RE_MAX = AU.STRICT_RE_MAX
P8_Q_UL_MIN = AU.STRICT_Q_UL_MIN''')
edit(B, '''                P8=dict(inlet_Re_limit=P8_RE_MAX, at_auditor_flow=p8(P8_Q_UL_MIN), at_design_flow=p8(p.Q_in_uL_min),
                        status='diagnostic only, not an audit and not gated (docs/LIMITATIONS.md 6): the two criteria '
                               'of CIF_v2\\'s physics auditor check P8, computed from the design (CIF_v2\\'s auditor is '
                               'not run on the exported mask; opencif\\'s own drawn-mask checks are '
                               'assessment.status.drawn_mask_assessment) at the auditor\\'s default operating point (5 uL/min) '
                               'and at the design flow; its 0.05 Re D_h is the entrance length spec v5, the glossary and '
                               'CIF_v2 S16-c call wrong in this regime, and its Re limit 1.0 is CIF_v2\\'s settable '
                               'default (App. B prints no value), not the design envelope\\'s (Re_max = 10, M14)'),''',
     '''                P8=dict(inlet_Re_limit=P8_RE_MAX, at_auditor_flow=p8(P8_Q_UL_MIN), at_design_flow=p8(p.Q_in_uL_min),
                        status='design-stage values (planned widths) of the physics auditor\\'s own inlet criteria, '
                               'check P8; the audit evaluates them on the drawn mask (validation.physics_audit, P8) '
                               'and gates the design envelope there (check M14: Re <= 10 on the channel, Re_gap and Re_delta '
                               '<= 1e-2); P8 is reported, not gated (docs/LIMITATIONS.md 6): its 0.05 Re D_h is the '
                               'entrance length spec v5, the glossary and CIF_v2 S16-c call wrong in this regime, and '
                               'its Re limit 1.0 is CIF_v2\\'s settable default (App. B prints no value), stricter than '
                               'the design envelope\\'s (Re_max = 10, M14)'),''')
edit(B, '''def build(spec, kind='linear', n_legs=1, outdir='.', name='device', failsafe=False, do_audit=True, verbose=True,
          R_turn=None, cut_tol=0.01, mode='strict', seal='require'):''',
     '''def build(spec, kind='linear', n_legs=1, outdir='.', name='device', failsafe=False, do_audit=True, verbose=True,
          R_turn=None, cut_tol=0.01, mode='strict', seal='require', ar_max=AR_MAX_DEFAULT):''')
edit(B, '''    rmin, where, nlt = min_curvature_radius(fluid)
    V['rounding'] = dict(min_local_radius_um=rmin, at=where, vertices_with_radius_below=nlt,
                         criterion_um=spec.rmin_ok, ok=bool(rmin >= spec.rmin_ok))
    say(''', '''    V['rounding'] = rounding_check(fluid, spec)
    rmin = V['rounding']['min_local_radius_um']
    say(''')
edit(B, '''    V['outlet_measured'] = {k: v for k, v in m['outlet'].items() if k != 'detail'}
    det_ = m['outlet'].get('detail') or {}
    # rc5: the leads' narrowest measured widths (the aspect-ratio inventory) and each side's stretch-end term
    V['outlet_measured']['lead_width_min_um'] = {k_: float(det_[k_]['w_lead_min']) for k_ in ('A', 'B', 'C')
                                                 if isinstance(det_.get(k_), dict) and 'w_lead_min' in det_[k_]}
    V['outlet_measured']['stretch_end'] = {k_: dict(c=float(det_[k_].get('c_end', 0.0)), R=float(det_[k_].get('R_end', 0.0)),
                                                    stretch_width_um=float(det_[k_]['stretch_ws']))
                                           for k_ in ('A', 'B') if isinstance(det_.get(k_), dict) and 'stretch_ws' in det_[k_]}
    V['loops_measured'] = [dict(row=int(lp['k']) + 1, A=dict(w_min=lp['A']['w_min'], R_rel=lp['A']['R'] / lp['target'] - 1),
                                B=dict(w_min=lp['B']['w_min'], R_rel=lp['B']['R'] / lp['target'] - 1),
                                Rc_rel_to_design=lp['Rc'] / lp['Rc_design'] - 1) for lp in m['loops']]
    if m['loops'] and 'side_continuity' in V:
        lwm = min(min(lp['A']['w_min'], lp['B']['w_min']) for lp in m['loops'])
        sc_v = V['side_continuity']
        sc_v['loop_narrowest_um'] = float(lwm)
        sc_v['narrowest_um'] = float(min(sc_v['narrowest_um'], lwm))
        sc_v['above_floor'] = bool(sc_v['narrowest_um'] >= spec.F_min - 0.05)
''', '''    outlet_and_loops_measured(V, m, spec)
''')
edit(B, '''    inten = intended_network(sc_f, info, m)
    f_star_rows = np.array([S.f_of(spec.delta_star, w) for w in sc_f.wc])
    turns_ = sc_f.info.get('turns', [])
    V['drawing_fidelity'] = drawing_fidelity(nw, inten, f_star_rows, sc_f.N, D_star)
    paper_, wR_, kL_ = L.transition_params(sc_f)
    nets = dict(full=inten)
    opts_ = active_options(spec, sc_f.f, f_star_rows, turns_, mode, seal)
    if OPT_PAPER in opts_:
        nets['no_paper'] = intended_network(sc_f, info, m, drop=('paper',))
    if OPT_PULL in opts_:
        nets['no_pull'] = intended_network(sc_f, info, m, drop=('pulled',))
    if OPT_PAPER in opts_ and OPT_PULL in opts_:
        nets['none'] = intended_network(sc_f, info, m, drop=('paper', 'pulled'))
    elif 'no_paper' in nets or 'no_pull' in nets:
        nets['none'] = nets.get('no_paper', nets.get('no_pull'))
    decl = declared_rows(spec, sc_f.f, f_star_rows, turns_, D_star, cut_tol, nets, mode, seal)
''', '''    inten, nets, decl = intended_and_declared(sc_f, info, m, spec, D_star, cut_tol, mode, seal)
    f_star_rows = np.array([S.f_of(spec.delta_star, w) for w in sc_f.wc])
    V['drawing_fidelity'] = drawing_fidelity(nw, inten, f_star_rows, sc_f.N, D_star)
    paper_, wR_, kL_ = L.transition_params(sc_f)
''')
edit(B, '''    first_lin = [a for k, s_, a, b in m['lin'] if k == kL and s_ == 'A']
    last_lin = [a for k, s_, a, b in m['lin'] if k == sc_f.N - 2 and s_ == 'A']
    wc_list = sorted(m['wc'])
    # the quantities the papers print, read back from the DRAWN mask (generated values, for comparison with the
    # printed ones in docs/REPRODUCTIONS.md; rc3 called this block 'published')
    V['drawn_parameters'] = dict(depth_um=spec.h, G_meas=float(np.median([a for k, s_, a, b in m['gaps'] if k >= sc_f.NM])),
                          G_meas_min=float(np.min([a for k, s_, a, b in m['gaps'] if k >= sc_f.NM])),
                          wc_in_meas=float(wc_list[0][1]), wc_out_meas=float(wc_list[-1][1]),
                          ws_first_linear_eq=float(first_lin[0]) if first_lin else None,
                          ws_out_eq=float(last_lin[0]) if last_lin else None,
                          footprint_mm=[(fluid.bounds[2] - fluid.bounds[0]) / 1000, (fluid.bounds[3] - fluid.bounds[1]) / 1000],
                          array_length_mm=float(sc_f.x[-1] / 1000),
                          head_length_mm=float((sc_f.x[sc_f.NM] - info['x_top']) / 1000) if sc_f.NM else 0.0)
''', '''    V['drawn_parameters'] = drawn_parameters(m, sc_f, fluid, info, kL)
''')
edit(B, '''    if do_audit:
        d_um = spec.F_min - 0.1
        aj = os.path.join(outdir, name + '_audit.json')
        if os.path.exists(aj):
            os.remove(aj)                                 # never read a previous run's result
        r = subprocess.run([sys.executable, DRC, aj, dxf, '%.2f' % d_um, '--islands'], capture_output=True, text=True)
        try:
            if r.returncode != 0:
                raise RuntimeError('DRC exited %d' % r.returncode)
            V['drc'] = json.load(open(aj))
        except Exception as e:
            V['drc'] = dict(error='%s: %s' % (e, r.stderr[-2000:]))
        say('[%s] DRC %s' % (name, {k: V['drc'].get(k) for k in ('n_components', 'invalid_rings', 'n_lost_open',
                                                                  'n_lost_solid', 'd_um')}))
''', '''    if do_audit:
        V['drc'] = run_drc(dxf, outdir, name, spec)
        say('[%s] DRC %s' % (name, {k: V['drc'].get(k) for k in ('n_components', 'invalid_rings', 'n_lost_open',
                                                                  'n_lost_solid', 'd_um')}))
    # 1.0.0: the physics auditor (P0-P10) on the photomask file read back: P0 requires it to equal the drawn polygon
    # vertex for vertex, so the measurement above is the measurement of the file
    V['physics_audit'] = AU.audit(dxf, sc_f, info, m, V, decl, mode, ar_max, fluid=fluid, nw=nw, intended=inten)
    say('[%s] %s (%.0fs)' % (name, AU.summary_line(V['physics_audit']), time.time() - t0))
''')

# the corner-radius check, as a function (the standalone audit evaluates it on the polygon read back)
edit(B, '''def components_drawn(fluid, sc, info):''', '''def rounding_check(fluid, spec):
    """V2: the smallest local radius of curvature of the drawn boundary against the spec's criterion"""
    rmin, where, nlt = min_curvature_radius(fluid)
    return dict(min_local_radius_um=rmin, at=where, vertices_with_radius_below=nlt, criterion_um=spec.rmin_ok,
                ok=bool(rmin >= spec.rmin_ok))


def outlet_and_loops_measured(V, m, spec):
    """the outlet and U-turn blocks of the record, from the measurement (and the loops' narrowest width into the
    side-channel continuity block)"""
    V['outlet_measured'] = {k: v for k, v in m['outlet'].items() if k != 'detail'}
    det_ = m['outlet'].get('detail') or {}
    # rc5: the leads' narrowest measured widths (the aspect-ratio inventory) and each side's stretch-end term
    V['outlet_measured']['lead_width_min_um'] = {k_: float(det_[k_]['w_lead_min']) for k_ in ('A', 'B', 'C')
                                                 if isinstance(det_.get(k_), dict) and 'w_lead_min' in det_[k_]}
    V['outlet_measured']['stretch_end'] = {k_: dict(c=float(det_[k_].get('c_end', 0.0)), R=float(det_[k_].get('R_end', 0.0)),
                                                    stretch_width_um=float(det_[k_]['stretch_ws']))
                                           for k_ in ('A', 'B') if isinstance(det_.get(k_), dict) and 'stretch_ws' in det_[k_]}
    V['loops_measured'] = [dict(row=int(lp['k']) + 1, A=dict(w_min=lp['A']['w_min'], R_rel=lp['A']['R'] / lp['target'] - 1),
                                B=dict(w_min=lp['B']['w_min'], R_rel=lp['B']['R'] / lp['target'] - 1),
                                Rc_rel_to_design=lp['Rc'] / lp['Rc_design'] - 1) for lp in m['loops']]
    if m['loops'] and 'side_continuity' in V:
        lwm = min(min(lp['A']['w_min'], lp['B']['w_min']) for lp in m['loops'])
        sc_v = V['side_continuity']
        sc_v['loop_narrowest_um'] = float(lwm)
        sc_v['narrowest_um'] = float(min(sc_v['narrowest_um'], lwm))
        sc_v['above_floor'] = bool(sc_v['narrowest_um'] >= spec.F_min - 0.05)


def drawn_parameters(m, sc_f, fluid, info, kL):
    """V8: the quantities the papers print, read back from the DRAWN mask (generated values, for comparison with the
    printed ones in docs/REPRODUCTIONS.md; rc3 called this block 'published')"""
    spec = sc_f.spec
    first_lin = [a for k, s_, a, b in m['lin'] if k == kL and s_ == 'A']
    last_lin = [a for k, s_, a, b in m['lin'] if k == sc_f.N - 2 and s_ == 'A']
    wc_list = sorted(m['wc'])
    return dict(depth_um=spec.h, G_meas=float(np.median([a for k, s_, a, b in m['gaps'] if k >= sc_f.NM])),
                G_meas_min=float(np.min([a for k, s_, a, b in m['gaps'] if k >= sc_f.NM])),
                wc_in_meas=float(wc_list[0][1]), wc_out_meas=float(wc_list[-1][1]),
                ws_first_linear_eq=float(first_lin[0]) if first_lin else None,
                ws_out_eq=float(last_lin[0]) if last_lin else None,
                footprint_mm=[(fluid.bounds[2] - fluid.bounds[0]) / 1000, (fluid.bounds[3] - fluid.bounds[1]) / 1000],
                array_length_mm=float(sc_f.x[-1] / 1000),
                head_length_mm=float((sc_f.x[sc_f.NM] - info['x_top']) / 1000) if sc_f.NM else 0.0)


def intended_and_declared(sc_f, info, m, spec, D_star, cut_tol, mode, seal):
    """the intended network (the closure plus the declared deliberate deviations) and the declared rows, both from the
    design and the drawing code's own profiles, before the drawn polygon is looked at.  Returns (inten, nets, decl)."""
    inten = intended_network(sc_f, info, m)
    f_star_rows = np.array([S.f_of(spec.delta_star, w) for w in sc_f.wc])
    turns_ = sc_f.info.get('turns', [])
    nets = dict(full=inten)
    opts_ = active_options(spec, sc_f.f, f_star_rows, turns_, mode, seal)
    if OPT_PAPER in opts_:
        nets['no_paper'] = intended_network(sc_f, info, m, drop=('paper',))
    if OPT_PULL in opts_:
        nets['no_pull'] = intended_network(sc_f, info, m, drop=('pulled',))
    if OPT_PAPER in opts_ and OPT_PULL in opts_:
        nets['none'] = intended_network(sc_f, info, m, drop=('paper', 'pulled'))
    elif 'no_paper' in nets or 'no_pull' in nets:
        nets['none'] = nets.get('no_paper', nets.get('no_pull'))
    decl = declared_rows(spec, sc_f.f, f_star_rows, turns_, D_star, cut_tol, nets, mode, seal)
    return inten, nets, decl


def run_drc(dxf, outdir, name, spec):
    """the DRC (opencif/checks/drc.py) on the photomask file, in its own process, at 0.1 um under the process floor"""
    d_um = spec.F_min - 0.1
    aj = os.path.join(outdir, name + '_audit.json')
    if os.path.exists(aj):
        os.remove(aj)                                 # never read a previous run's result
    r = subprocess.run([sys.executable, DRC, aj, dxf, '%.2f' % d_um, '--islands'], capture_output=True, text=True)
    try:
        if r.returncode != 0:
            raise RuntimeError('DRC exited %d' % r.returncode)
        return json.load(open(aj))
    except Exception as e:
        return dict(error='%s: %s' % (e, r.stderr[-2000:]))


def components_drawn(fluid, sc, info):''')
# ---------------------------------------------------------------------------------------------------------- pipeline.py
P = 'opencif/pipeline.py'
edit(P, '''                    length_budget_G9=('NOT_GIVEN', 'WITHIN', 'OVER'), drawn_mask_assessment=('RUN',),
                    physics_auditor_P8=('NOT_RUN_ON_MASK',), fabricated_or_measured=('NO',),''',
     '''                    length_budget_G9=('NOT_GIVEN', 'WITHIN', 'OVER'), drawn_mask_assessment=('RUN',),
                    physics_audit=('PASS', 'PASS_OUTSIDE_ENVELOPE', 'FAIL', 'INCOMPLETE'), fabricated_or_measured=('NO',),''')
edit(P, '''    if rc['declared_rows']:
        chk('declared departures (reported, not a pass condition)', True,''',
     '''    # 1.0.0: the physics auditor on the photomask file (opencif.checks.auditor): its hard checks are a pass condition in
    # every mode; its P8 (the Reynolds and Stokes gates on the drawn mask) is part of the validity envelope
    pa = V.get('physics_audit')
    if pa is None:
        chk('physics audit of the photomask file (P0-P10)', False, 'not run', 'the audit must run', 'incomplete')
    else:
        drc_missing = not V.get('drc') or 'error' in V.get('drc', {})
        chk('physics audit of the photomask file (P0-P10)', not pa['hard_failures'],
            dict(status=pa['status'], failed=pa['hard_failures'],
                 checks={c['id']: ('ok' if c['ok'] else ('reported' if c['role'] == 'report' else 'FAIL'))
                         for c in pa['checks']}),
            'the photomask read back from its file passes the physics auditor\\'s hard checks P0, R0, P1-P6b, P9, R6b '
            'and P10, and P9-AR in strict mode (CIF_v2\\'s P1-P10 ported; P0: the file is the drawing; R0: the '
            'auditor\\'s own ruler on the straight rows; R6b: the posts; P7 and P8 are reported, M14 is part of the '
            'validity envelope; validation.physics_audit has each check, its value and its requirement)',
            'incomplete' if (drc_missing and pa['hard_failures'] == ['P9']) else 'hard')
        p8_ = next(c for c in pa['checks'] if c['id'] == 'M14')
        chk('physics audit M14: the design code\\'s Reynolds and Stokes gates on the drawn mask', p8_['ok'],
            dict((k, p8_['value'][k]) for k in ('Re_channel_max', 'Re_gap_max', 'Re_delta_max')), p8_['requirement'],
            'envelope')
    if rc['declared_rows']:
        chk('declared departures (reported, not a pass condition)', True,''')
edit(P, '''    env_txt = ('OUTSIDE THE VALIDITY ENVELOPE (%s): the cutoff law is not calibrated there; %s'
               % ('; '.join(why), band_txt(band)))''', '''    if pa is not None and env['Reynolds'] and env['Stokes'] and not p8_['ok']:
        why.append('on the drawn mask (physics audit M14): Re %.3g, Re_gap %.2e, Re_delta %.2e (limits %g, %g)'
                   % (p8_['value']['Re_channel_max'], p8_['value']['Re_gap_max'], p8_['value']['Re_delta_max'],
                      BU.RE_MAX, BU.RE_STOKES))
    env_txt = ('OUTSIDE THE VALIDITY ENVELOPE (%s): the cutoff law is not calibrated there; %s'
               % ('; '.join(why), band_txt(band)))''')
edit(P, '''                  physics_auditor_P8='NOT RUN ON THE MASK: CIF_v2\\'s physics auditor (P1-P10) is not run on the exported '
                                     'mask; its P8 criteria are computed from the design as diagnostics, not gated '
                                     '(validation.gates.P8)',''',
     '''                  physics_audit=_audit_status(V.get('physics_audit'), V.get('drc')),''')
edit(P, '''                 drawn_mask_assessment='RUN', physics_auditor_P8='NOT_RUN_ON_MASK', fabricated_or_measured='NO',''',
     '''                 drawn_mask_assessment='RUN', physics_audit=_audit_code(V.get('physics_audit'), V.get('drc')),
                 fabricated_or_measured='NO',''')
edit(P, '''def assess(rep, spec, opts, fluid, notes, chip_rect):''', '''def _audit_code(pa, drc):
    if pa is None:
        return 'INCOMPLETE'
    if pa['hard_failures'] == ['P9'] and (not drc or 'error' in drc):
        return 'INCOMPLETE'
    return pa['status']


def _audit_status(pa, drc):
    """the physics audit's status line (assessment.status.physics_audit)"""
    code = _audit_code(pa, drc)
    if pa is None:
        return 'INCOMPLETE: the physics audit did not run'
    p8s = next((c for c in pa['checks'] if c['id'] == 'P8'), None)
    strict_txt = (' The auditor\\'s own, stricter inlet limit (P8, Re <= 1; reported): %s.' % p8s['value']['result']) if p8s else ''
    if code == 'INCOMPLETE':
        return 'INCOMPLETE: the DRC did not run, so P9 could not be evaluated' + strict_txt
    if code == 'FAIL':
        return 'FAIL: %s (validation.physics_audit)' % ', '.join(pa['hard_failures']) + strict_txt
    return ('%s: the photomask file read back passes the physics auditor (CIF_v2\\'s P1-P10 ported; P0: the file is '
            'the drawing; R0: its straight rows measured again by the auditor\\'s own ruler; R6b: its posts are the '
            'process\\'s)%s.' % (code, ' with M14 outside the validity envelope' if code == 'PASS_OUTSIDE_ENVELOPE'
                                 else '') + strict_txt)


def assess(rep, spec, opts, fluid, notes, chip_rect):''')
# the aspect-ratio inventory: one function, shared with the auditor's P9
edit(P, '''    dp = V.get('drawn_parameters') or {}
    widths = {'gap opening G': float(dp.get('G_meas_min', spec.G))}
    if dp.get('wc_in_meas') and dp.get('wc_out_meas'):
        widths['centre channel'] = min(float(dp['wc_in_meas']), float(dp['wc_out_meas']))
    if sc_:
        widths['straight side channel'] = sc_['narrowest_um']
    if lw:
        widths['U-turn side channel'] = min(lw)
    if cd is not None and cd.get('head_unit_widths_um'):
        widths['meander head unit'] = min(cd['head_unit_widths_um'])
    leads = list(((V.get('outlet_measured') or {}).get('lead_width_min_um') or {}).values())
    if leads:
        widths['outlet lead'] = min(leads)
''', '''    widths = AU.drawn_widths(V, spec)
''')
edit(P, "from . import build as BU", "from . import build as BU\nfrom .checks import auditor as AU")
edit(P, "    lw = [min(lp['A']['w_min'], lp['B']['w_min']) for lp in V.get('loops_measured', [])]\n", '')
edit(P, '''        rep, fluid, sc_f, info = BU.build(spec, kind, n_legs, outdir, name, failsafe=opts['failsafe'],
                                          do_audit=opts['drc'], verbose=verbose, R_turn=opts['leg_spacing'],
                                          cut_tol=opts['cut_tol'], mode=opts['mode'], seal=opts['seal'])''',
     '''        rep, fluid, sc_f, info = BU.build(spec, kind, n_legs, outdir, name, failsafe=opts['failsafe'],
                                          do_audit=opts['drc'], verbose=verbose, R_turn=opts['leg_spacing'],
                                          cut_tol=opts['cut_tol'], mode=opts['mode'], seal=opts['seal'],
                                          ar_max=opts.get('AR_max', BU.AR_MAX_DEFAULT))''')
print('patch_audit: %d edits' % n_done)
