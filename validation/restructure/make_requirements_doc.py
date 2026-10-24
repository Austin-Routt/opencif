"""Write docs/CIF_V2_REQUIREMENTS.md for opencif 1.0.0 from the read-only inventory of CIF_v2.py's requirements
(compiled on the 1.0.0 working tree before the design code's reports were added). The rows that those reports, R6b,
and the new tests changed are updated, the counts are recomputed from the rows, the prose is written in the author's
voice, and the closing sections are rewritten.
usage: python -I make_requirements_doc.py INVENTORY_MD OUT_MD"""
import re
import sys

src, out = sys.argv[1], sys.argv[2]
s = open(src, encoding='utf-8').read()

DCR = '`validation.design_code_reports`'
ROWS = {   # row id prefix -> (status, level, evidence)
    '| B3 total pressure drop ': ('PORTED-REPORTED', 'build',
        'design_code.py `reports` B3: %s.B3_pressure_drop: trunk_kPa (the design code\'s dP_tot = sum R_c(i) Q_c(i), '
        'its upstream-flow convention) and drawn_device_kPa (the inlet pressure of the drawn network, gap 1 to the '
        'outlets), at the design flow and the spec\'s fluid' % DCR),
    '| B6 F_real = eta F ': ('PORTED-REPORTED', 'build',
        '%s.B6_F_real: eta F_drawn with CIF_v2\'s eta 0.903, reported, never applied (LIMITATIONS 5); trace: eta = '
        'CIF_v2 Pairing.eta' % DCR),
    '| B7 transfer bias (S5, P3) ': ('PORTED-REPORTED', 'build',
        '%s.B7_transfer_bias: [min, max] sqrt(mean c2(beta_cal)/c2(h/w_c(i))) over the rows, never applied; '
        'kernels.c2, c2_exact and Pairing.transfer_bias bit-identical to CIF_v2\'s (test_design_code.py)' % DCR),
    '| B8 sensitivity map (S16-e) ': ('PORTED-REPORTED', 'build',
        '%s.B8_sensitivity_map: design_code.sensitivity_map, CIF_v2\'s finite differences (d = 1e-4) at the outlet '
        'and gap 1, = CIF_v2.sensitivity_map within 1e-6 on its reference device (test_design_code.py)' % DCR),
    '| B9 wall porosity phi = G/Lambda ': ('PORTED-REPORTED', 'build', '%s.B9_wall_porosity: phi = G/(G + g_p)' % DCR),
    '| P3 transfer bias reported, never applied ': ('PORTED-REPORTED', 'build',
        'as B7: the band is reported on every build, D is never corrected'),
    '| M9 invariance sweeps ': ('PORTED-GATED', 'test',
        'test_design_code.py: the schedule (f, P, S, widths, closure, N) identical at 3, 30 and 300 uL/min and '
        '1-2.5 mPa s, and N, F, f, P, S, w_c identical at post lengths 20 and 40 um (the pitch column); the GS column '
        'concerns CIF_v2\'s folded side channels, which opencif does not draw'),
    '| M17 sensitivity map emitted ': ('PORTED-REPORTED', 'build',
        '%s.B8_sensitivity_map.M17_widths_carry_the_device (None when an outlet channel is not narrower than h); '
        'reported, as every design_code report' % DCR),
    '| R1 determinism ': ('PORTED-GATED', 'test',
        'test_design_code.py: one specification drawn in two separate processes and in the test\'s own gives the '
        'same photomask, byte for byte; the release rebuild compares every photomask with RC8\'s (docs/VALIDATION.md)'),
    '| R6b DXF posts = the 2:1 pill ': ('PORTED-GATED', 'build',
        'auditor R6b (hard): every post of every straight row read back from the photomask file, across within '
        '0.5 um of l_g, along within 0.5 um of g_p, ratio within 0.05 (CIF_v2\'s three criteria, on every post '
        'instead of one; seeded fault in test_auditor.py)'),
    '| M19 Gifford Fig. 2, engine-driven ': ('PORTED-GATED', 'test',
        'test_design_code.py, section 6: from opencif\'s planner at F = 2/3 (the endpoint w_s = w_c), N = ceil(i_T) and '
        'f(i) = f_gap for the six f_gap = k x 9.6e-5 of panel A; Q_s/Q_c equal to [(1-2f)^-i - 1]/2 within 1e-12; at '
        'f_gap = 1.76e-4 the six depths give one N and one flow-ratio curve while w_s fans out; the six k x 9.6e-5 end '
        'within 2% of the curves digitized from the published figure. CIF_v2.gifford_fig2\'s own f_list (5.76, 4.83, '
        '3.90, 2.97, 2.04, 0.96) x 1e-4 is not the caption\'s linear sequence and ends 5% early on the fifth curve '
        '(an erratum of the frozen design code; its M19 checks internal consistency only); the tutorial draws the four '
        'panels (section 4.5)'),
    '| Z3 embedded run_selftest passes ': ('EQUIVALENT', 'test',
        'test_auditor.py, part 3: fourteen seeded faults on a drawn device, caught by P4 (an intended network 1% off), '
        'P0, P1, P4 and P5 (a widened gap), P3, R6b, R0, P4 and P5 (a stepped side wall), P0 with R0 and R6b, M14, P9, '
        'P10, P6, and P9-AR, each the check meant to catch its fault'),
    '| V_op Re ': ('PORTED-REPORTED', 'build',
        '%s.V_operating_regime at 5 uL/min and at the design flow (= CIF_v2.V_operating_regime on its reference '
        'device: test_design_code.py); gates.P8 (planned widths) and auditor P8 (the drawn inlet)' % DCR),
    '| V_op entrance ': ('PORTED-REPORTED', 'build', 'as V_op Re'),
    '| V_op shear ': ('PORTED-REPORTED', 'build',
        '%s.V_operating_regime.*.gdot_per_s = C(alpha) Q/(a b^2) at the inlet, against 1500 1/s, at 5 uL/min and '
        'the design flow (gates.shear_peak: every row at the design flow)' % DCR),
    '| V_op pressure ': ('PORTED-REPORTED', 'build',
        '%s.V_operating_regime.*.dP_kPa: the design code\'s trunk-only sum r(w_c) Lambda P Q against 200 kPa; B3 '
        'for the drawn device' % DCR),
    '| V_steric D/G ': ('PORTED-REPORTED', 'build',
        '%s.V_steric.D_over_G against 0.5 (= CIF_v2.V_steric: test_design_code.py); passage_report: the measured gap '
        'opening against 2 D (`validation.passage.classes`)' % DCR),
    '| V_steric D/h ': ('PORTED-REPORTED', 'build', '%s.V_steric.D_over_h against 0.3' % DCR),
    '| V_steric delta*/h ': ('PORTED-REPORTED', 'build',
        '%s.V_steric.dstar_over_h against 0.2; stricter and gated as E11 delta*/h <= 0.10 (A9)' % DCR),
    '| V_fab_tolerance ': ('PORTED-REPORTED', 'build',
        '%s.V_fab_tolerance: 24 draws of N(0, 0.25 um) on every centre and side width (seed 0, the design code\'s '
        'order), the design network solved again; yield = share with |F - F_target| <= 3e-3, ok iff >= 0.95, '
        'adherence mean/min; deterministic and exact at sigma 0 (test_design_code.py); a report, as there' % DCR),
}
lines = s.split('\n')
done = set()
for i, ln in enumerate(lines):
    for pre, (st, lev, ev) in ROWS.items():
        if ln.startswith(pre):
            cells = ln.split(' | ')
            # | id | lines | what | role | status | level | evidence |
            assert len(cells) >= 7, ln
            cells[4] = '**%s**' % st
            cells[5] = lev
            lines[i] = ' | '.join(cells[:6] + [ev + ' |'])
            done.add(pre)
missing = set(ROWS) - done
assert not missing, missing
s = '\n'.join(lines)

# recount: the rows of sections A, B, D, E, F, G
sec_of, counts = None, {}
STAT = ['PORTED-GATED', 'PORTED-REPORTED', 'EQUIVALENT', 'NOT-APPLICABLE', 'MISSING']
for ln in s.split('\n'):
    m = re.match(r'^## ([A-G])\. ', ln)
    if m:
        sec_of = m.group(1)
        continue
    if ln.startswith('## '):
        sec_of = None
        continue
    if sec_of in ('A', 'B', 'D', 'E', 'F', 'G') and ln.startswith('| ') and not ln.startswith('| ID') \
            and not ln.startswith('|---'):
        m = re.search(r'\*\*(PORTED-GATED|PORTED-REPORTED|EQUIVALENT|NOT-APPLICABLE|MISSING)\*\*', ln)
        if m:
            counts.setdefault(sec_of, dict.fromkeys(STAT, 0))[m.group(1)] += 1
tab = ['| section | rows | PORTED-GATED | PORTED-REPORTED | EQUIVALENT | NOT-APPLICABLE | MISSING |',
       '|---|---|---|---|---|---|---|']
tot = dict.fromkeys(STAT, 0)
for sec in ('A', 'B', 'D', 'E', 'F', 'G'):
    c = counts[sec]
    tab.append('| %s | %d | %s |' % (sec, sum(c.values()), ' | '.join(str(c[k]) for k in STAT)))
    for k in STAT:
        tot[k] += c[k]
tab.append('| **all** | **%d** | %s |' % (sum(tot.values()), ' | '.join('**%d**' % tot[k] for k in STAT)))
a = s.index('| section | rows |')
b = s.index('\n\n', a)
s = s[:a] + '\n'.join(tab) + s[b:]

# the opening paragraph
a = s.index('Read-only audit, 2026-10-22.')
b = s.index('\n\n', a)
s = s[:a] + (
    'This document lists every gate, refusal, mandatory report, validator, and acceptance criterion of the design code '
    '`reference/CIF_v2.py` v2.1.0 (8,058 lines, md5 0e0961a5..., asserted by tests/test_kernels.py) and of its embedded '
    '`_PHYSICS_SRC` (ladder_physics.py) and `_MEANDER_SRC` (meander_pattern.py), read in full, with the status of each '
    'in opencif 1.0.0. CIF_v2 was also imported read-only to list its registry, in which `VALIDATORS` holds 162 '
    'entries. The inventory was compiled on 2026-10-22 on the 1.0.0 working tree. The design code\'s reports (section B '
    'and the V_* validators), R6b, and the tests of M9, M19, and R1 were added to 1.0.0 afterward, and their rows give '
    'their final status. opencif references are given as file:line in code {{STAMP}}.') + s[b:]
s = s.replace('"?" marks a status I am not certain of (listed again at the end).',
              '"?" marks a status that is a judgment (listed again at the end).')
assert '"?" marks a status that is a judgment' in s

# the closing sections
a = s.index('## Every MISSING item')
s = s[:a] + """## What 1.0.0 added

The inventory of the working tree found eleven requirements with no counterpart in opencif. Every one of them is a
report or an acceptance-only validator in CIF_v2, and none is a pass condition there. They are the five mandatory
reports of `design()` that opencif lacked (B3 pressure drop, B6 F_real, B7 transfer bias, B8 sensitivity map, and B9
porosity), the validators that test two of them (P3 and M17), and four parts of the acceptance run's operating
validators (the shear and pressure limits of `V_operating_regime`, D/h of `V_steric`, and `V_fab_tolerance`). opencif
1.0.0 computes each by the design code's own formula on every build and records it without gating it, as the design
code does, in `validation.design_code_reports` (`opencif/checks/design_code.py`, with each formula checked against
CIF_v2's own function by `tests/test_design_code.py`). It also turned four rows that the inventory could not mark with
certainty into ported checks. R6b (the DXF's posts are the process's pill) is the auditor's hard check R6b on every
post of the straight rows of the file, and M9 (invariance), M19 (Gifford's Fig. 2 driven by the engine), and R1
(determinism) are tests (`tests/test_design_code.py`). No row is MISSING.

## Item whose status is a judgment

G24, the depth-correct field solve (CIF_v2, line 3606), is marked EQUIVALENT. opencif has no two-dimensional field
solve on the mask. A short across a wall is caught by the hole count (a cut post changes it), by the per-gap
measurement of the width, by the auditor's R0 and R6b, and by Q1.

## Documentation corrected in 1.0.0

The inventory also found statements in opencif's documentation that the 1.0.0 working tree no longer supported. They
were the statement that the auditor is not run on the exported mask (LIMITATIONS 6, ALIGNMENT AB27, DISSERTATION_NOTES,
EVALUATION, the reference README, and foundation/README's D95), the statement that its other checks have only
counterparts of purpose, the statement that `V_steric` and `V_fab_tolerance` are not computed, the module names of RC8
(`cifsched` and `measure2`), and an O(f^2) argument in a kernel's docstring that ALIGNMENT AG9 says opencif does not
use. Each is corrected in 1.0.0 (`docs/CHANGELOG.md`, 1.0.0).

## Method notes

Each row was checked against the source text of both trees, and where opencif's counterpart is a test, the level
column says so. A check of the same formula and threshold that is a pass condition or a refusal in opencif is
PORTED-GATED even where CIF_v2 only tests it on its own fixtures, because its validate() battery never sees a user's
design.
"""

# The prose blocks of the inventory, in the author's voice. The rows of the tables keep their terse form.
PROSE = [
    ("""## How CIF_v2 checks things (this decides what "port" can mean)

* **On the design being built** CIF_v2 enforces only section A: `design()` and `march()` raise on the plateau, gamma_GS, w_min, the lamina and taper limits, M13 and M14; the constructors and mask builders raise on bowtie posts, impossible U-turns, routability and (v2) composition interference; the unified compiler adds four pre-flight gates (section E, CU1-CU4).  `design()` also emits mandatory reports (section B).
* **`validate()`** (162 validators, section D) is a self-test of the implementation on CIF_v2's own fixtures (`_Ctx`, `reference_request()`, demo specs A/B/C).  It never sees a user's design.  Its own printout says so: "This proves the IMPLEMENTATION is correct.  It cannot prove the MODEL is."
* **`run_acceptance()`** passes when `goldens and G14 and jsuite` (8041-8043) hold on CIF_v2's optimizer devices; the D-suite, `V_operating_regime`, `V_steric` and `V_fab_tolerance` are written to its JSON but are not in `ok`.  `V_outlet_divider`, `J5`, `J10` and `self_test()` are defined and never called.
* **The physics auditor** (`audit()`, P1-P10 with P6b) runs on CIF_v2's meander demo patterns; Z1, Z2, T1, W3, W4, Q1, Q5 and R2a require every P check to pass, and `compile_unified` refuses when P3 or P4 fails.""",
     """## How CIF_v2 applies its checks

How CIF_v2 applies a check decides what a port of it can mean, and CIF_v2 applies its checks in four ways. First, on
the design being built, CIF_v2 enforces only section A. `design()` and `march()` raise on the plateau, gamma_GS, w_min,
the lamina and taper limits, M13, and M14; the constructors and mask builders raise on bowtie posts, impossible U-turns,
routability, and (in v2) composition interference; and the unified compiler adds four pre-flight gates (section E,
CU1 to CU4). `design()` also emits mandatory reports (section B). Second, `validate()` (162 validators, section D) is a
self-test of the implementation on CIF_v2's own fixtures (`_Ctx`, `reference_request()`, and the demo specifications
A, B, and C), and it never sees a user's design, as its own printout states ("This proves the IMPLEMENTATION is
correct.  It cannot prove the MODEL is."). Third, `run_acceptance()` passes when `goldens and G14 and jsuite` (lines
8041 to 8043) hold on CIF_v2's optimizer devices, whereas the D-suite, `V_operating_regime`, `V_steric`, and
`V_fab_tolerance` are written to its JSON without entering `ok`; `V_outlet_divider`, `J5`, `J10`, and `self_test()` are
defined and never called. Finally, the physics auditor (`audit()`, P1 to P10 with P6b) runs on CIF_v2's meander demo
patterns, the validators Z1, Z2, T1, W3, W4, Q1, Q5, and R2a require every P check to pass, and `compile_unified`
refuses when P3 or P4 fails."""),
    ("""Statuses: **PORTED-GATED** (the same check is a pass condition or refusal in opencif), **PORTED-REPORTED** (the same check is computed and recorded, never gated), **EQUIVALENT** (same purpose, different method: the evidence says how), **NOT-APPLICABLE** (concerns geometry or machinery opencif does not have, or CIF_v2's own fixtures), **MISSING**.  "level": where opencif's counterpart lives: build (every design), test (tests/), both, or -.  "?" marks a status that is a judgment (listed again at the end).""",
     """Each row carries one of five statuses. PORTED-GATED means that the same check is a pass condition or a refusal in
opencif, PORTED-REPORTED that the same check is computed and recorded without being gated, EQUIVALENT that the check
serves the same purpose by a different method (the evidence says how), NOT-APPLICABLE that it concerns geometry or
machinery that opencif does not have, or CIF_v2's own fixtures, and MISSING that opencif has no counterpart. The column
"level" says where opencif's counterpart lives, which is the build (every design), the tests (tests/), both, or
neither (-). A "?" marks a status that is a judgment, and such rows are listed again at the end."""),
    ("""CIF_v2 roles: REFUSE (raises on the request), TRIP (AssertionError tripwire), REPORT (Device.reports/report()), VAL (validate() pass condition on CIF_v2's fixtures), PRE (compile_unified pre-flight refusal), AUDIT (audit() check), ACC-PASS (in run_acceptance's ok), ACC-REP (acceptance JSON only), DEF (defined, never called), SELF (embedded self-test).""",
     """The column of CIF_v2's roles uses ten labels, which are REFUSE (raises on the request), TRIP (an AssertionError
tripwire), REPORT (Device.reports or report()), VAL (a pass condition of validate() on CIF_v2's fixtures), PRE (a
pre-flight refusal of compile_unified), AUDIT (a check of audit()), ACC-PASS (in run_acceptance's ok), ACC-REP (in the
acceptance JSON only), DEF (defined and never called), and SELF (the embedded self-test)."""),
    ("""Section A+B are what CIF_v2 enforces or reports on a design: 49 rows.  The 47 supplied constants (section C) are tallied separately below, since most of them are the thresholds of the rows above.""",
     """Sections A and B, with 49 rows, hold what CIF_v2 enforces or reports on a design. The 47 supplied constants
(section C) are tallied separately below, because most of them are the thresholds of the rows above."""),
    ("""These are the only CIF_v2 checks that act on a user's request in `design()` / `build_mask()`: everything in section D runs on CIF_v2's own fixtures.""",
     """These are the only CIF_v2 checks that act on a user's request in `design()` and `build_mask()`, because everything in
section D runs on CIF_v2's own fixtures."""),
    ("""The trace test accounts for every Tolerances/Pairing/Process field (test_appendix_trace.py:330-379; ALIGNMENT sec. 5).  It does not account for the validators, reports and acceptance functions of the other sections.""",
     """The trace test accounts for every field of Tolerances, Pairing, and Process (test_appendix_trace.py:330-379; ALIGNMENT,
section 5). It does not account for the validators, reports, and acceptance functions of the other sections."""),
    ("""Constants by disposition: consumed 3, gated 8, not used / not applicable 11, partly gated 1, ported to the auditor 5, reported 1, spec input, replaced or reused by another check 18.""",
     """By disposition, 3 constants are consumed, 8 are gated, 11 are not used or not applicable, 1 is partly gated, 5 are
ported to the auditor, 1 is reported, and 18 are specification inputs, replaced, or reused by another check."""),
    ("""Role "VAL" = a pass condition of `validate()`, which runs on CIF_v2's own fixtures (`_Ctx`: a short FOLD device, its DIRECT twin, `reference_request()`, the embedded demo specs), never on a user's design.  An opencif counterpart in the build is therefore stronger than the original; a counterpart in `tests/` is the like-for-like.""",
     """The role "VAL" is a pass condition of `validate()`, which runs on CIF_v2's own fixtures (`_Ctx`, which is a short FOLD
device and its DIRECT twin, `reference_request()`, and the embedded demo specifications) and never on a user's design.
An opencif counterpart in the build is therefore stronger than the original, whereas a counterpart in `tests/` is the
like-for-like."""),
    ("""run_acceptance's pass condition is `ok = goldens and G14 and jsuite` (8041-8043); the D-suite, V_operating_regime, V_steric and V_fab_tolerance are stored in the acceptance JSON but are not part of `ok`.""",
     """The pass condition of run_acceptance is `ok = goldens and G14 and jsuite` (lines 8041 to 8043), whereas the D-suite,
V_operating_regime, V_steric, and V_fab_tolerance are stored in the acceptance JSON without being part of `ok`."""),
    ("""## F. The embedded physics core (_PHYSICS_SRC = ladder_physics.py, CIF_v2 line 3718; line numbers below are in the extracted copy aud/_PHYSICS_SRC.py)""",
     """## F. The embedded physics core (_PHYSICS_SRC, which is ladder_physics.py, at CIF_v2 line 3718; the line numbers below are those of the module extracted from that string)"""),
    ("""Role "AUDIT" = a check of `audit()`; validators Z1, Z2, T1, W3, W4, Q1, Q5 and R2a require every P check to pass on the demo devices, and compile_unified refuses when P3 or P4 fails (U3).  opencif's port is opencif/checks/auditor.py, run on every photomask (build.py:1332-1335) and made a pass condition by pipeline.assess (pipeline.py:679-697).""",
     """The role "AUDIT" is a check of `audit()`. The validators Z1, Z2, T1, W3, W4, Q1, Q5, and R2a require every P check to
pass on the demo devices, and compile_unified refuses when P3 or P4 fails (U3). opencif's port is
opencif/checks/auditor.py, which runs on every photomask (build.py:1332-1335) and is made a pass condition by
pipeline.assess (pipeline.py:679-697)."""),
    ("""## G. The embedded meander core (_MEANDER_SRC = meander_pattern.py, CIF_v2 line 3716)""",
     """## G. The embedded meander core (_MEANDER_SRC, which is meander_pattern.py, at CIF_v2 line 3716)"""),
    ("""Its 19 check_* functions are the Y1-Y19 rows of section D (one validator per check); they are not repeated here.""",
     """Its 19 check_* functions are the rows Y1 to Y19 of section D (one validator per check), so they are not repeated here."""),
]
for old, new in PROSE:
    assert s.count(old) == 1, old[:80]
    s = s.replace(old, new)
s = s.replace('P9-posts', 'R6b').replace('auditor P8-M14', 'auditor M14').replace('P8-M14', 'M14')
s = s.replace('# CIF_v2.py requirements and their status in opencif 1.0.0',
              '# The requirements of CIF_v2.py and their status in opencif 1.0.0')
s = re.sub(r'\*\*([^*\n]+)\*\*', r'\1', s)                           # no bold, in the tables either
s = s.replace(' -- ', '; ').replace('judgement', 'judgment')
# American spelling in the prose of the rows; code spans and identifiers keep their own
AMERICAN = {'centre': 'center', 'centres': 'centers', 'centreline': 'centerline', 'centrelines': 'centerlines',
            'neighbour': 'neighbor', 'neighbours': 'neighbors', 'neighbouring': 'neighboring',
            'equalised': 'equalized', 'equalises': 'equalizes', 'minimised': 'minimized', 'minimises': 'minimizes',
            'parameterisation': 'parameterization', 'realised': 'realized', 'realises': 'realizes',
            'realise': 'realize', 'modelled': 'modeled', 'behaviour': 'behavior', 'behaviours': 'behaviors',
            'normalised': 'normalized', 'labelled': 'labeled', 'micrometres': 'micrometers', 'grey': 'gray',
            'optimised': 'optimized', 'recognised': 'recognized', 'characterised': 'characterized',
            'generalised': 'generalized', 'generalises': 'generalizes', 'summarised': 'summarized',
            'analysed': 'analyzed', 'colour': 'color', 'colours': 'colors'}


def _american(seg):
    seg = re.sub(r'(\d) %(?![\w%])', r'\1%', seg)                     # a closed-up percent sign
    return re.sub(r'(?<![\w`])(%s)(?![\w`])' % '|'.join(sorted(AMERICAN, key=len, reverse=True)),
                  lambda m: AMERICAN[m.group(1)], seg)


s = ''.join(part if i % 2 else _american(part) for i, part in enumerate(re.split(r'(`[^`\n]*`)', s)))
# The line references to opencif's code and tests, resolved against the final tree (inventory/LINEREF_MAP.json).
import json as _json                                                   # noqa: E402
import os as _os                                                       # noqa: E402
_map_path = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'inventory', 'LINEREF_MAP.json')
n_ref = 0
if _os.path.exists(_map_path):
    _map = _json.load(open(_map_path, encoding='utf-8'))

    def _sub(mo):
        global n_ref
        new = _map.get(mo.group(0))
        if new is not None and new != mo.group(0):
            n_ref += 1
            return new
        return mo.group(0)
    s = re.sub(r'(?<![\w./])[a-z_0-9]+\.py:\d+(?:-\d+)?(?:,\d+(?:-\d+)?)*', _sub, s)
# Three rows whose evidence the map cannot correct, because a reference lacks its file name or the row's text no longer
# describes the code at the reference (the line-reference review of the final tree).
AFTER_MAP = [
    ('pipeline.py:805 env `foot`, 1249-1250 env_bad', 'pipeline.py:805 env `foot`, pipeline.py:1333-1334 env_bad'),
    ('schedule.py:804-805 (400000 rows) -> SpecError (pipeline.py:1051-1062)',
     'schedule.py:804-805 raises RuntimeError after 400000 rows, as CIF_v2 does, and pipeline.prepare converts it into '
     'a SpecError (pipeline.py:1051-1062), so the command line reports a spec error (exit 6); seeded in '
     'test_negative.py:177-187'),
    ('auditor KCL residual < 1e-9 (auditor.py:546)',
     'auditor P2 KCL residual below 1e-4 of the smallest undeclared gap flow (auditor.py:531-546)'),
]
if _os.path.exists(_map_path):
    for _old, _new in AFTER_MAP:
        if s.count(_old) != 1:
            raise SystemExit('make_requirements_doc: %r found %d times' % (_old, s.count(_old)))
        s = s.replace(_old, _new)
open(out, 'w', encoding='utf-8').write(s)
print('wrote %s: %s; %d line references moved' % (out, ', '.join('%s %d' % (k, v) for k, v in tot.items()), n_ref))
