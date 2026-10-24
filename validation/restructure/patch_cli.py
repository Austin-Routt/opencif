"""opencif 1.0.0: the `audit` command, the teaching core and the new test suites.  Applied after patch_audit.py.
usage: python -I patch_cli.py TREE TOOLS_DIR"""
import os
import shutil
import sys

T, HERE = sys.argv[1:3]
n_done = 0


def edit(rel, old, new, count=1):
    global n_done
    p = os.path.join(T, rel)
    s = open(p, encoding='utf-8').read()
    assert s.count(old) == count, (rel, s.count(old), old[:100])
    s = s.replace(old, new)
    open(p, 'w', encoding='utf-8').write(s)
    n_done += 1


MN = 'opencif/__main__.py'
edit(MN, """   python -m opencif check SPEC.yaml          (validate the spec, plan the design and evaluate the gates; no drawing;
                                              exit 5 refused, 3 outside the envelope, 0 inside, 6 spec error)
""", """   python -m opencif check SPEC.yaml          (validate the spec, plan the design and evaluate the gates; no drawing;
                                              exit 5 refused, 3 outside the envelope, 0 inside, 6 spec error)
   python -m opencif audit RUN_DIR [--name NAME] [--out FILE] [--no-drc]
                                             (the physics auditor on a finished build's photomask file, from scratch:
                                              plan and draw again, read the DXF back, measure it, check P0-P10,
                                              R0, M14 and R6b; exit 0 pass, 3 pass with M14 outside the envelope,
                                              2 fail, 4 incomplete)
""")
edit(MN, """    c.add_argument('spec', help='a spec YAML or JSON file')
""", """    c.add_argument('spec', help='a spec YAML or JSON file')
    au = sub.add_parser('audit', help='run the physics auditor (P0-P10, R0, M14, R6b) on a finished build\\'s photomask '
                                      'file, from scratch (exit 0 pass, 3 pass with M14 outside the envelope, 2 fail, '
                                      '4 incomplete)')
    au.add_argument('run_dir', help='the build\\'s folder (its record <name>.json and photomask <name>.dxf)')
    au.add_argument('--name', help='the build\\'s name, when the folder holds more than one record')
    au.add_argument('--out', help='also write the audit as JSON to this file')
    au.add_argument('--no-drc', action='store_true', help='skip the DRC (P9 is then incomplete, exit 4)')
""")
edit(MN, """def _main(a):
    if a.cmd == 'generate':""", """def _audit(a):
    \"\"\"opencif audit RUN_DIR: the physics auditor on a finished build, from its photomask file\"\"\"
    from .checks import auditor as AU
    r = AU.audit_run(a.run_dir, a.name, drc=not a.no_drc)
    for c in r['checks']:
        print('%-6s %-9s %-8s %s' % (c['id'], 'ok' if c['ok'] else ('reported' if c['role'] == 'report' else 'FAIL'),
                                     c['role'], c['check']))
    if a.out:
        with open(a.out, 'w') as fh:
            json.dump(r, fh, indent=1, default=float)
        print('audit: %s' % a.out)
    status = 'INCOMPLETE' if (a.no_drc and r['hard_failures'] == ['P9']) else r['status']
    print(AU.summary_line(r) + ('' if status == r['status'] else ' (INCOMPLETE: the DRC was not run)'))
    return {'PASS': 0, 'PASS_OUTSIDE_ENVELOPE': 3, 'FAIL': 2, 'INCOMPLETE': 4}[status]


def _main(a):
    if a.cmd == 'audit':
        return _audit(a)
    if a.cmd == 'generate':""")
edit('tests/test_rc7.py', "    'drawn_mask_assessment', 'physics_auditor_P8', 'fabricated_or_measured', 'reproduction', 'source_trace'})",
     "    'drawn_mask_assessment', 'physics_audit', 'fabricated_or_measured', 'reproduction', 'source_trace'})")
CO = 'opencif/compose.py'
edit(CO, """Checks on the drawn array: DRC at the floor; one fluid component (shared inlet) or n (separate); in every element the""",
     """The element must be a clean PASS whose photomask passed the physics auditor (its record's validation.physics_audit;
opencif.checks.auditor): the copies are its geometry, translated.
Checks on the drawn array: DRC at the floor; one fluid component (shared inlet) or n (separate); in every element the""")
edit(CO, """    bm = (rec.get('design', {}).get('spec', {}) or {}).get('bend_model', 'stokes')""",
     """    pa = (rec.get('validation') or {}).get('physics_audit') or {}
    if pa.get('status') != 'PASS':
        return False, ('only an element whose photomask passed the physics auditor can be arrayed (validation.physics_audit '
                       'status %r)' % pa.get('status'))
    bm = (rec.get('design', {}).get('spec', {}) or {}).get('bend_model', 'stokes')""")
edit(CO, """    out = dict(name=aname, element=name, element_verdict=verdict, element_code=el_code, code=code, count=count, pitch_um=pitch,""",
     """    out = dict(name=aname, element=name, element_verdict=verdict, element_code=el_code,
               element_physics_audit=(rec.get('validation') or {}).get('physics_audit', {}).get('status'),
               code=code, count=count, pitch_um=pitch,""")
edit('tests/test_units.py', """# 8. compose takes only a clean pass
for v, c, want in (('PASS', 0, True), ('REPRODUCED with declared departures (seed)', 0, False),
                   ('PASS -- OUTSIDE THE VALIDITY ENVELOPE (u_min 0.03)', 3, False), ('PASS with an allowed departure (seed)', 0, False),
                   ('FAIL (Q1)', 2, False)):
    check('compose element gate: %s' % v[:40], CO.element_ok(dict(assessment=dict(verdict=v, exit_code=c)))[0] == want)""",
     """# 8. compose takes only a clean pass whose photomask passed the physics auditor
for v, c, pa, want in (('PASS', 0, 'PASS', True), ('REPRODUCED with declared departures (seed)', 0, 'PASS', False),
                       ('PASS -- OUTSIDE THE VALIDITY ENVELOPE (u_min 0.03)', 3, 'PASS', False),
                       ('PASS with an allowed departure (seed)', 0, 'PASS', False), ('FAIL (Q1)', 2, 'FAIL', False),
                       ('PASS', 0, None, False), ('PASS', 0, 'PASS_OUTSIDE_ENVELOPE', False)):
    rec_ = dict(assessment=dict(verdict=v, exit_code=c), validation=dict(physics_audit=dict(status=pa)) if pa else {})
    check('compose element gate: %s, audit %s' % (v[:40], pa), CO.element_ok(rec_)[0] == want)""")
edit('opencif/api.py', """    @property
    def ok(self):
        \"\"\"exit code 0: every hard check of the generator passed inside the validity envelope -- a statement about the
        1-D model on the drawn mask, not that the chip is fabricable or that a device was made or measured.""",
     """    @property
    def ok(self):
        \"\"\"exit code 0: every hard check of the generator passed inside the validity envelope (the physics audit of the
        photomask file among them) -- a statement about the 1-D model on the drawn mask, not that the chip is fabricable
        or that a device was made or measured.""")
edit('opencif/api.py', """    @property
    def over_chip_area(self):""", """    @property
    def physics_audit(self):
        \"\"\"the physics audit of the photomask file (opencif.checks.auditor): 'PASS', 'PASS_OUTSIDE_ENVELOPE' or
        'FAIL', or None for a refused design (no photomask); every check is in record['validation']['physics_audit']\"\"\"
        return (((self.record or {}).get('validation') or {}).get('physics_audit') or {}).get('status')

    @property
    def over_chip_area(self):""")
edit('opencif/pipeline.py', """    if a.get('status_codes'):                              # rc7 (ADJ-C6)""",
     """    pa = V.get('physics_audit')
    if pa:                                                 # 1.0.0: the physics audit, check by check
        lines += ['', '**Physics audit of %s: %s** (opencif.checks.auditor; each check\\'s value and requirement are in the '
                  'record, validation.physics_audit)' % (pa['photomask']['file'], pa['status']), '',
                  '| check | | result | role |', '|---|---|---|---|']
        for c in pa['checks']:
            lines.append('| %s | %s | %s | %s |' % (c['id'], c['check'], 'ok' if c['ok'] else
                                                    ('reported' if c['role'] == 'report' else '**fail**'), c['role']))
    if a.get('status_codes'):                              # rc7 (ADJ-C6)""")
edit('opencif/pipeline.py', """        drc_['file_final'] = dxf_name
    rep['files'] = dict(""", """        drc_['file_final'] = dxf_name
    pa_ = (rep.get('validation') or {}).get('physics_audit')
    if isinstance(pa_, dict) and (pa_.get('photomask') or {}).get('file') not in (None, dxf_name):
        pa_['photomask']['file_final'] = dxf_name         # 1.0.0: the audited file (the same bytes), under the verdict's name
    rep['files'] = dict(""")
edit('opencif/pipeline.py', """            drc_['file_final'] = fn_                      # rc8: as run()
""", """            drc_['file_final'] = fn_                      # rc8: as run()
        pa_ = (rep.get('validation') or {}).get('physics_audit')
        if isinstance(pa_, dict) and (pa_.get('photomask') or {}).get('file') not in (None, fn_):
            pa_['photomask']['file_final'] = fn_          # 1.0.0: as run()
""")
edit('opencif/pipeline.py', """'record, validation.physics_audit)' % (pa['photomask']['file'], pa['status']), '',""",
     """'record, validation.physics_audit)' % (pa['photomask'].get('file_final', pa['photomask']['file']), pa['status']), '',""")
# the design code's reports and acceptance validators (opencif/checks/design_code.py), recorded in every build
shutil.copy(os.path.join(HERE, 'design_code.py'), os.path.join(T, 'opencif', 'checks', 'design_code.py'))
# the physics auditor's own ruler (opencif/checks/ruler.py): the straight rows measured again from the file alone
shutil.copy(os.path.join(HERE, 'ruler.py'), os.path.join(T, 'opencif', 'checks', 'ruler.py'))
edit('opencif/build.py', "from opencif.checks import auditor as AU  # noqa: E402\n",
     "from opencif.checks import auditor as AU  # noqa: E402\nfrom opencif.checks import design_code as DCR  # noqa: E402\n")
edit('opencif/build.py', """    say('[%s] %s (%.0fs)' % (name, AU.summary_line(V['physics_audit']), time.time() - t0))
""", """    say('[%s] %s (%.0fs)' % (name, AU.summary_line(V['physics_audit']), time.time() - t0))
    # 1.0.0: the design code's mandatory reports and acceptance validators, each by its own formula (reported, as there)
    V['design_code_reports'] = DCR.reports(sc_f, info, m, spec, float(nw['F']))
""")
edit('opencif/pipeline.py', """    if rc['declared_rows']:
        chk('declared departures (reported, not a pass condition)', True,""",
     """    dcr = V.get('design_code_reports')
    if dcr:                                               # 1.0.0: reported, as in the design code
        vo_ = dcr['V_operating_regime']['at_5_uL_min']
        chk('the design code\\'s reports and acceptance validators (reported, not a pass condition)', True,
            dict(dP_device_kPa=dcr['B3_pressure_drop']['drawn_device_kPa'], eta_F=dcr['B6_F_real']['eta_F'],
                 transfer_bias=dcr['B7_transfer_bias']['D_true_over_D_predicted'],
                 V_operating_regime_at_5_uL_min=vo_['ok'], V_steric=dcr['V_steric']['ok'],
                 V_fab_tolerance_yield=dcr['V_fab_tolerance']['yield_F']),
            'design()\\'s mandatory reports (pressure drop, eta F, transfer bias, sensitivity map, wall porosity) and '
            'run_acceptance\\'s operating validators (V_operating_regime at 5 uL/min, V_steric, V_fab_tolerance), each by '
            'the design code\\'s formula: validation.design_code_reports', 'model')
    if rc['declared_rows']:
        chk('declared departures (reported, not a pass condition)', True,""")
# 1.0.0 final: a plan that cannot be completed (for example, a march that needs more than 400,000 rows) is a spec
# error (exit 6) rather than a traceback (exit 1), and the figure title uses a colon where RC8 used an em dash.
edit('opencif/pipeline.py', """        return _prepare(spec, opts, say)
    except SP.SpecError:
        raise
    except ValueError as e:""", """        return _prepare(spec, opts, say)
    except SP.SpecError:
        raise
    except (ValueError, RuntimeError) as e:""")
edit('opencif/pipeline.py', """        title = '%s \u2014 D %.3f µm, F %.4f;""", """        title = '%s: D %.3f µm, F %.4f;""")
shutil.copy(os.path.join(HERE, 'teaching.py'), os.path.join(T, 'opencif', 'teaching.py'))
for t in ('test_kernels.py', 'test_auditor.py', 'test_teaching.py', 'test_design_code.py'):
    shutil.copy(os.path.join(HERE, t), os.path.join(T, 'tests', t))
print('patch_cli: %d edits, teaching.py and 3 test suites copied' % n_done)
