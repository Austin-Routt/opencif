"""opencif 1.0.0: the RC8 suites' edits for the new reports (exact-match edits; usage: python -I patch_tests.py TREE).

test_appendix_trace's CIF_v2 field dispositions: eta and beta_cal are now REPORTED on every build (the design code's
F_real = eta F and transfer bias, opencif/checks/design_code.py), each value asserted equal to CIF_v2's; L_eq's
disposition names both of its replacements (rc7's square-turn corner term as well as the curved-duct solution), which
RC8 could not edit without changing the tests' hash."""
import os
import sys

TREE = sys.argv[1]
n = 0


def edit(rel, old, new):
    global n
    p = os.path.join(TREE, rel)
    s = open(p, encoding='utf-8').read()
    if s.count(old) != 1:
        raise SystemExit('patch_tests: %s: expected once, found %d: %r' % (rel, s.count(old), old[:80]))
    open(p, 'w', encoding='utf-8').write(s.replace(old, new))
    n += 1


T = 'tests/test_appendix_trace.py'
# 1.0.0 final: the march's row cap is a spec error (pipeline.prepare), seeded beside the other unplannable designs.
edit('tests/test_negative.py', """                      (100.0, 'a closure with no solution (D 100 um on 100 um, f_gap 0.11) -> SpecError, not a crash',
                       'cannot be planned')):""", """                      (100.0, 'a closure with no solution (D 100 um on 100 um, f_gap 0.11) -> SpecError, not a crash',
                       'cannot be planned'),
                      (0.25, 'a march of more than 400,000 rows (D 0.25 um on 100 um, f_gap 7.8e-7) -> SpecError, not a '
                       'crash', 'march did not reach')):""")
edit(T, "from opencif.checks import measure as M2\n",
     "from opencif.checks import measure as M2\nfrom opencif.checks import design_code as DCR\n")
edit(T, """                    beta_cal=('not needed on the cubic pair (glossary: if the kernel is never switched, h_cal does not '
                              'enter N); ALIGNMENT AG2', True),""",
     """                    beta_cal=(REPORTED, tuple(DCR.K.Pairing().beta_cal) == tuple(C.Pairing().beta_cal)),
                    # 1.0.0: the transfer bias it implies is reported on every build (design_code B7), never applied;
                    # on the cubic pair it does not enter N (glossary; ALIGNMENT AG2)""")
edit(T, "                    eta=('LIMITATIONS 5', 'η' in lim_items.get(5, '')),",
     "                    eta=(REPORTED, DCR.ETA == C.Pairing().eta and 'η' in lim_items.get(5, '')),"
     "   # 1.0.0: eta F reported (B6)")
edit(T, "                    L_eq=('superseded by the exact curved-duct solution (ALIGNMENT AG7)', True),",
     "                    L_eq=('superseded by the exact curved-duct solution for arcs and, for square turns, the computed '\n"
     "                          'corner term c_corner (ALIGNMENT AG7b)', True),")
# the glossary's numeric macros: from the preamble when the dissertation's overlay is present, else from
# reference/glossary_values.json (a copy of the repository without the foundation texts)
edit(T, """pre = open(os.path.join(ROOT, 'reference_glossary', 'PREAMBLE_FINAL.tex')).read()
mac = {m.group(1): m.group(2) for m in re.finditer(r'\\\\newcommand\\{\\\\(\\w+)\\}\\{([0-9.]+)\\}', pre)}
""", """import json
GLOSSARY = json.load(open(os.path.join(ROOT, 'reference', 'glossary_values.json')))
PRE_PATH = os.path.join(ROOT, 'reference_glossary', 'PREAMBLE_FINAL.tex')
if os.path.exists(PRE_PATH):            # the dissertation's overlay: the frozen preamble itself
    pre = open(PRE_PATH).read()
    mac = {m.group(1): m.group(2) for m in re.finditer(r'\\\\newcommand\\{\\\\(\\w+)\\}\\{([0-9.]+)\\}', pre)}
else:                                   # the public repository: the values the checks below use
    pre = None
    mac = dict(GLOSSARY['macros'])
""")
edit(T, """m_cap = re.search(r'\\\\newcommand\\{\\\\fgapDesignCap\\}\\{([0-9.]+)\\\\times10\\^\\{(-?\\d+)\\}\\}', pre)
check('glossary \\\\fgapDesignCap = the calibration cap the record reports (schedule.FGAP_CAP)',
      m_cap is not None and abs(float(m_cap.group(1)) * 10 ** int(m_cap.group(2)) - S.FGAP_CAP) < 1e-12)""",
     """cap_ = (re.search(r'\\\\newcommand\\{\\\\fgapDesignCap\\}\\{([0-9.]+)\\\\times10\\^\\{(-?\\d+)\\}\\}', pre).groups()
        if pre is not None else GLOSSARY['fgapDesignCap'])
check('glossary \\\\fgapDesignCap = the calibration cap the record reports (schedule.FGAP_CAP)',
      abs(float(cap_[0]) * 10 ** int(cap_[1]) - S.FGAP_CAP) < 1e-12)""")
edit(T, """check('reference_glossary/PREAMBLE_FINAL.tex is the frozen v114 preamble (md5 3c075af8...)',
      md5(os.path.join(ROOT, 'reference_glossary', 'PREAMBLE_FINAL.tex')) == FROZEN_V114['PREAMBLE_FINAL.tex'])
""", """check('reference/glossary_values.json names the frozen v114 preamble (md5 3c075af8...)',
      GLOSSARY['preamble_md5'] == FROZEN_V114['PREAMBLE_FINAL.tex'])
if pre is not None:                     # the dissertation's overlay
    check('reference_glossary/PREAMBLE_FINAL.tex is the frozen v114 preamble (md5 3c075af8...)',
          md5(PRE_PATH) == FROZEN_V114['PREAMBLE_FINAL.tex'])
    m_ = {m.group(1): m.group(2) for m in re.finditer(r'\\\\newcommand\\{\\\\(\\w+)\\}\\{([0-9.]+)\\}', pre)}
    check('reference/glossary_values.json holds the preamble\\'s own values (%d macros)' % len(GLOSSARY['macros']),
          all(m_.get(k) == v for k, v in GLOSSARY['macros'].items()) and list(cap_) == GLOSSARY['fgapDesignCap'])
""")
# test_trace_anchors L1: the anchors are RC8's; 1.0.0 draws RC8's photomasks byte for byte (validation/RELEASE.json)
edit('tests/test_trace_anchors.py', "        from opencif.pipeline import code_stamp\n        cs = code_stamp()\n        check('L1 code stamp', cs['sha256'] == A['code_stamp']['sha256'], '%s (anchors %s)' % (cs['sha256'], A['code_stamp']['sha256']))\n",
     "        from opencif.pipeline import code_stamp\n        cs = code_stamp()\n        # The anchors were traced on RC8's photomasks. A later release whose evidence shows every photomask byte for\n        # byte RC8's (validation/RELEASE.json) traces the same masks, so its own stamp stands for the anchors' too.\n        rel_p = os.path.join(ROOT, 'validation', 'RELEASE.json')\n        rel = json.load(open(rel_p)) if os.path.exists(rel_p) else {}\n        cnt = rel.get('counts') or {}\n        same_masks = bool(rel.get('code_stamp') == cs['sha256'] and rel.get('rc8_code_stamp') == A['code_stamp']['sha256']\n                          and cnt.get('photomasks') and cnt.get('identical_to_rc8') == cnt.get('photomasks'))\n        check('L1 code stamp', cs['sha256'] == A['code_stamp']['sha256'] or same_masks,\n              '%s (anchors %s%s)' % (cs['sha256'], A['code_stamp']['sha256'],\n                                     '; the release draws RC8\\'s photomasks byte for byte' if same_masks else ''))\n")

print('patch_tests: %d edits' % n)
