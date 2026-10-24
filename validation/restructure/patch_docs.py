"""opencif 1.0.0: the documents' edits for the release (exact-match edits; usage: python -I patch_docs.py TREE TOOLS_V1).

What changes: the physics auditor now runs on every photomask with its own ruler (LIMITATIONS 6 and 12, ALIGNMENT AB27,
DISSERTATION_NOTES § 4, foundation/README's D95, the errata list's D95 and D106); the design code's reports are
recorded (LIMITATIONS 5, ALIGNMENT AG7i and § 5, § 4's eta); RC8's module names in the documents that describe the
current code (the history documents keep the names their releases had: CHANGELOG, REVIEW, provenance); the RC8 README
becomes docs/REFERENCE.md beside a new README; the new documents are copied in."""
import os
import re
import shutil
import sys

TREE, HERE = sys.argv[1], sys.argv[2]
n = 0


def edit(rel, old, new, count=1):
    global n
    p = os.path.join(TREE, rel)
    s = open(p, encoding='utf-8').read()
    if s.count(old) != count:
        raise SystemExit('patch_docs: %s: expected %d, found %d: %r' % (rel, count, s.count(old), old[:90]))
    open(p, 'w', encoding='utf-8').write(s.replace(old, new))
    n += 1


def between(rel, start, end, new):
    """replace the text from `start` (included) to `end` (excluded); each must occur once"""
    global n
    p = os.path.join(TREE, rel)
    s = open(p, encoding='utf-8').read()
    if s.count(start) != 1 or s.count(end) != 1:
        raise SystemExit('patch_docs: %s: markers %d, %d' % (rel, s.count(start), s.count(end)))
    a, b = s.index(start), s.index(end)
    open(p, 'w', encoding='utf-8').write(s[:a] + new + s[b:])
    n += 1


# ---------------------------------------------------------------------------------------------------------------------
# docs/LIMITATIONS.md
# ---------------------------------------------------------------------------------------------------------------------
LIM = 'docs/LIMITATIONS.md'
edit(LIM, """5. **F is the 1-D filtrate fraction.**  3-D flow at the outlet bifurcation (the glossary's η = 0.903, a measured
   correction) is not applied; Dinh's network F is 10.7 % above his measurement by construction (App. A R8; App. B
   P-C4.2).
""", """5. **F is the 1-D filtrate fraction.**  3-D flow at the outlet bifurcation (the glossary's η = 0.903, a measured
   correction) is not applied; Dinh's network F is 10.7 % above his measurement by construction (App. A R8; App. B
   P-C4.2).  Since 1.0.0 every record reports ηF beside F, as the design code does (*record:*
   `validation.design_code_reports.B6_F_real`; reported, never applied).
""")
edit(LIM, """6. **CIF_v2's E11 and M16 are gates (rc3); the physics auditor P8 is reported, not gated.** rc3 ports CIF_v2's""",
     """6. **CIF_v2's E11 and M16 are gates (rc3); the physics auditor runs on every photomask (1.0.0), and its stricter P8
   is reported, not gated.** rc3 ports CIF_v2's""")
between(LIM, "   *The physics auditor.*  App. B gives", "7. **Aspect ratio (rc4: a gate in strict mode).**",
        """   *The physics auditor (1.0.0: ported, and run on every photomask).*  App. B gives the design code the channel gate
   Re_max = 10 (M14) and says that "the physics auditor separately checks the drawn device's inlet against its own,
   stricter channel limit" (P8). CIF_v2 has both: `design()` gates Re ≤ 10, and its physics path requires inlet
   Re ≤ 1.0 (a settable default; App. B gives no value) at the device's flow (`Q_in_uLmin`, default 5 µL/min): the
   auditor's check `P8_stokes_gates`, which adds a naive entrance length 0.05 Re D_h ≤ Λ/2, and the unified compiler's
   pre-flight Reynolds gate; its validators Z1–Z3 and T1 require its demo devices, at 2 and 5 µL/min, to pass. Since
   1.0.0 opencif runs that auditor, ported (`opencif/checks/auditor.py`, `docs/PHYSICS_AUDIT.md`), on the photomask
   file of every build, read back from the DXF.  P0 (the file is the drawing), R0 (its own ruler on the straight rows,
   below), P1–P6b, P9, R6b (the posts) and P10 are pass conditions in every mode (P7 reported, as in CIF_v2); the
   design code's own gate M14 (Re ≤ 10 on the channel, Re_gap and Re_δ ≤ 10⁻²) is evaluated on the drawn mask as the
   check M14 (outside it the verdict is exit 3).  P8 itself is checked on the drawn inlet at the design flow and at
   5 µL/min and recorded under its own name (*record:* `validation.physics_audit`, check P8; `gates.P8` gives the
   planned values), but it is not a pass condition.  At their own flows (25–800 µL/min) every opencif design exceeds
   P8's Re limit (inlet Re 3.6–5.6 on the 27 built strict designs, 6.7 on Gifford Fig. 6, 61 on Xia RS-7), so a P8
   gate would refuse them all, including the Dinh V1 presets at 30 µL/min (inlet Re 4.65), an operating point App. B's
   scope discussion puts inside the design code's channel gate (worst-case channel Re ≈ 5.2, against Re_max = 10);
   CIF_v2 runs its own lineage demo device at 2 µL/min. At 5 µL/min every design passes both criteria (inlet Re
   0.38–0.93; naive L_e at most 0.18 of Λ/2). And the L_e criterion uses the formula CIF_v2's own `entrance_length`
   calls "simply the wrong formula" in this regime, which spec v5 also rejects as a gate, and the glossary calls it
   "wrong here".  The port keeps every check's purpose and threshold; where a DXF file shows something other than
   CIF_v2's pattern objects, the formula is the one the file supports: P1 anchors the cubic pair opencif uses (CIF_v2's
   physics module anchors its laminar kernel), P6b uses App. A's march–outlet identity (CIF_v2's corrected closed form
   exists only for its canonical taper), and P9, R6b and P10 are measured on the polygon (a DXF has no strokes or
   fillet primitives).  `tests/test_auditor.py` runs the port and CIF_v2's own auditor on CIF_v2's three demonstration
   devices and requires the same numbers for the network, the terminal window, P3 (its PE part), P4, P5, P7 and P8 and
   the same verdicts, and seeds a fault for every check.  App. B's sentence that the C7/C8 software is "warranted by
   their acceptance suite and the physics auditor" (:4375–4377) refers back to "the physics auditor" of :4349–4350,
   which App. B's provenance tag at :4371–4372 lists beside `CIF_v2.py`'s `reynolds` and `design` (M14).  Up to RC8
   opencif did not run it (*record:* `assessment.status.physics_auditor_P8`, "NOT RUN ON THE MASK"), and open item D95
   proposed naming CIF_v2 there.  In 1.0.0 the sentence holds of opencif's auditor, with one difference the
   dissertation should state: the auditor's stricter limit is checked and reported, the design code's M14 is enforced.
   D95 is therefore closed into D106, revised for 1.0.0, which names opencif 1.0.0 and says which gate its auditor
   enforces (`docs/DISSERTATION_NOTES.md` § 4; the frozen text is not edited here).

   *The auditor's own ruler.*  The generator measures its drawing where its plan puts each part and takes the parts'
   lengths from the plan (`opencif.checks.measure`); an auditor using only that measurement would share its errors.
   The auditor therefore measures every ordinary straight row again (`opencif.checks.ruler`, which shares no code with
   the generator's measurement): the posts found among the file's islands, the gap openings, the wall thickness, the
   row pitch, the centre width at the post tips and the side width to the outer wall, all from the file; the two must
   agree within 10⁻³ (check R0), and the auditor's network uses the ruler's resistances there.  The head units, the
   U-turns and the outlet are not measured by the ruler: they stay on the generator's measurement, model terms and
   plan-taken lengths (item 12).

   *CIF_v2's mandatory reports and acceptance validators (1.0.0: reported on every build).*  CIF_v2's `design()`
   attaches reports to every device (the total pressure drop, F_real = ηF, the transfer bias, the sensitivity map, the
   wall porosity), and its acceptance run (`run_acceptance`) computes an operating-regime contract
   (`V_operating_regime`: Re ≤ 1.0 at 5 µL/min, 0.05 Re D_h ≤ Λ/2, wall shear ≤ 1500 s⁻¹, ΔP ≤ 200 kPa), a steric
   screen (`V_steric`: D/G < 0.5, D/h < 0.3, δ\\*/h < 0.2) and a Monte-Carlo width tolerance (`V_fab_tolerance`:
   σ = 0.25 µm on every centre and side width, 24 draws, yield = the share of draws with |F − F_target| ≤ 3 × 10⁻³, ok
   when ≥ 0.95); none of them is a pass condition there (`ok` = the gold devices, G14 and the J-suite). Since 1.0.0
   opencif computes each by the design code's own formula on every build and records it, also without making it a
   pass condition (*record:* `validation.design_code_reports`; `opencif/checks/design_code.py`;
   `tests/test_design_code.py` checks each formula against CIF_v2's own function).  On the shipped designs:
   {{DCR_LIM6}}  `docs/CIF_V2_REQUIREMENTS.md` lists every gate, validator and acceptance criterion of `CIF_v2.py` and
   its status in opencif.
""")
edit(LIM, """12. **Parts taken from the design, not measured on the drawing:** the post length across the flow (ℓ_g) and so the side
    channels' inner edge, the radii of head and outlet arcs (their widths are measured), the head rows' centre-path
    lengths, off-chip tubing.  *Record:* notes list them.
""", """12. **Parts taken from the design, not measured on the drawing:** the radii of head, U-turn and outlet arcs (their
    widths are measured), the head rows' centre-path lengths and the lengths of the generator's straight segments,
    off-chip tubing.  *Record:* notes list them.  Since 1.0.0 the physics audit measures the straight rows' lengths and
    the post length across the flow (ℓ_g, which every gap resistance takes from the design) on the photomask file
    itself: its own ruler reads the row pitch and the wall thickness from the posts and holds the generator's
    measurement to them (R0), and R6b requires every straight-row post's extent across the flow to be ℓ_g within
    0.5 µm.  The head units, the U-turns and the outlet are not re-measured by the auditor.
""")
edit(LIM, "the curved-duct solution's (`cifsched.depth_lambda`:", "the curved-duct solution's (`schedule.depth_lambda`:")
edit(LIM, "Since rc5 the closure counts it (`cifsched.outlet_exp_c`; item 8)", "Since rc5 the closure counts it (`schedule.outlet_exp_c`; item 8)")
edit(LIM, "(`ZeroDivisionError` in `measure2`)", "(`ZeroDivisionError` in `measure2`, now `opencif.checks.measure`)")
edit(LIM, "(`cifsched` lines 441, 548, 572,", "(`opencif/physics/schedule.py` lines 441, 548, 572,")
edit(LIM, "`cifsched` \"opencif F10\"", "`schedule` \"opencif F10\"")

# ---------------------------------------------------------------------------------------------------------------------
# foundation/README.md and the errata list (published with the repository, beside the frozen texts they describe)
# ---------------------------------------------------------------------------------------------------------------------
FR = 'foundation/README.md'
between(FR, "**Open item D95** (recorded 2026-10-10 in the rc4 review round", "**Open items D96 and D97**",
        """**Open item D95, closed into D106 (opencif 1.0.0).**  D95 (recorded 2026-10-10 in the rc4 review round, wording
revised in the rc5 round) concerned App. B :4373–4377: the targets-to-photomask tool and the machine-enforced envelope
(Prospectus C7, C8) "are delivered as verified software, warranted by their acceptance suite and the physics auditor
rather than by a defeasible inference".  "The physics auditor" refers back to :4349–4350 ("the physics auditor
separately checks the drawn device's inlet against its own, stricter channel limit"), whose item closes with the tag
"code: CIF\\_v2.py reynolds, design (M14); physics auditor (P8); reference implementation, VISCOSITY\\_PLASMA"
(:4371–4372).  The auditor is CIF_v2's (its `_PHYSICS_SRC` carries the P8 check).  Up to RC8 opencif did not run it on
its exported mask, so the sentence was true of CIF_v2 and not of opencif, and D95 proposed naming the program: "...
warranted by their acceptance suite (and, in the design code \\texttt{CIF\\_v2.py}, its physics auditor) ...".
opencif 1.0.0 runs the auditor, ported, on the photomask file of every build (`docs/PHYSICS_AUDIT.md`): its P8 checks
the drawn inlet against the auditor's own stricter limit and records the result, and the design code's channel and
Stokes gates (M14) are enforced on the drawn mask.  The sentence now holds of opencif's auditor as well, with that
difference, which the erratum states rather than leaves to the reader: D95 is closed into **D106, revised for 1.0.0**
(`ERRATA_NEXT_UNFREEZE.md`).  The note on "reference implementation" stands: in the frozen texts that phrase names the
preprint's notebook (App. B :4562; glossary bibliography, "CIF reference implementation v3.20", a computational
notebook; App. B :3492–3493 cites "reference implementation \\S 8.7, pe\\_consistency\\_check"), and the tag's
VISCOSITY\\_PLASMA is that program's constant (App. B :4359; `CIF_v2.py` has no VISCOSITY\\_PLASMA), so the
frozen texts' "reference implementation" is not `CIF_v2.py`, which App. A calls "the design code". The sentence
carries no argument node (App. B: "The silence is a decision"), so no claim, label, gate or number moves.

""")
ER = 'foundation/ERRATA_NEXT_UNFREEZE.md'
edit(ER, """## D95 \u2014 App. B :4373–4377, the warrant sentence (wording)

Recorded in the rc4 round, wording revised in rc5; the full text is in `README.md` ("Open item D95").  No number moves.
""", """## D95 \u2014 App. B :4373–4377, the warrant sentence (wording): closed into D106

Recorded in the rc4 round, wording revised in rc5; closed into D106 in opencif 1.0.0, which runs the physics auditor on
every photomask (`README.md`, "Open item D95").  No number moves.

## D106 \u2014 name the generator beside the design code (proposed with RC8; revised for opencif 1.0.0)

**Places.**  App. A's title page and table caption (A:50, A:184–187), which name the programs behind the results; App.
B :4373–4377, the software sentence, which speaks of "the targets-to-photomask tool (Prospectus C7)" without naming it.

**Reason.**  The dissertation's software is opencif; the foundation names only the design code `CIF_v2.py`.  opencif
1.0.0 restates the design code's kernels (its `opencif.physics.kernels`, tested to return `CIF_v2.py`'s numbers bit
for bit) and runs the design code's physics auditor, ported, on every photomask; the auditor enforces the design code's
channel and Stokes gates (M14) on the drawn mask and checks and reports, without enforcing, its own stricter inlet limit
(P8: inlet Re ≤ 1, which every lineage device exceeds at its own flow).

**Proposed text.**  A:184–187, after "the design code `CIF_v2.py` (Results 10–14)": "whose kernels the targets-to-
photomask generator opencif 1.0.0 (code stamp `{{STAMP}}`) reproduces bit for bit".  B:4375–4377: "... warranted by
their acceptance suite and the physics auditor (the design code's, which opencif 1.0.0 runs on every photomask file:
it enforces the design code's channel and Stokes gates on the drawn mask and reports its own stricter inlet limit)
rather than by a defeasible inference."  No claim, label, gate or number moves.
""")
# The texts are published with the repository, so that every citation of them resolves in it: the PDFs compiled at the
# freeze and the differences from the previous freeze (v110) that docs/ALIGNMENT.md read go beside the LaTeX sources
# that RC8's tree carries, each checked against its md5.
import hashlib                                                    # noqa: E402
FOUNDATION_ADD = {'APPENDIX_A_revised.pdf': 'a972b6e3b1c595cb4ffffd175f542d13',
                  'APPENDIX_B_current.pdf': 'af17d9db8481b17ba2fcfe3f70a70e4f',
                  'glossary_official.pdf': '9b656101ac6f00b14d2e57164a18e0b5',
                  'v110_to_v114_APPENDIX_A.diff': '49b6fa487cc52533e788590122c8365f',
                  'v110_to_v114_APPENDIX_B.diff': '1c4289ea0d8f17f802fdaacddb90f07c',
                  'v110_to_v114_glossary.diff': 'e610791c85e4b810aba88fa15c972162'}
for f, h in FOUNDATION_ADD.items():
    src = os.path.join(HERE, 'foundation_add', f)
    if hashlib.md5(open(src, 'rb').read()).hexdigest() != h:
        raise SystemExit('patch_docs: foundation_add/%s is not the frozen file (md5 %s...)' % (f, h[:8]))
    shutil.copy(src, os.path.join(TREE, 'foundation', f))
edit(FR, """# The frozen foundation (v114, 2026-09-29)

Copies of the sources opencif is aligned against (`docs/ALIGNMENT.md`).  The dissertation's own copies govern; these are
included so the evaluation archive is self-contained.  Appendices A and B and the glossary converged and were frozen at
v114 on 2026-09-29 (the project's freeze record, `CIF_ARGUMENTATION_FROZEN_v114.tar.gz`, md5
`c6abcca40fce4dcfdf38e5e6471a5a1a`).  `tests/test_appendix_trace.py` asserts the md5s of the four tex files below when
this folder is present.
""", """# The texts that opencif cites (frozen at v114, 2026-09-29)

This folder holds the author's texts that the code and the documents of opencif cite. Appendix A gives the step-by-step
derivations of the design equations (Results 1 to 15). Appendix B gives the argument for the design method, with its
thesis, its claims C1 to C9, their premises, and their critical questions. The glossary is the reference for every
variable, unit, and equation, and its preamble `PREAMBLE_FINAL.tex` holds the glossary's numeric values as macros. The
three texts were written as the mathematical foundation of the author's dissertation and were frozen on 2026-09-29
(v114; the freeze record `CIF_ARGUMENTATION_FROZEN_v114.tar.gz` has md5 `c6abcca40fce4dcfdf38e5e6471a5a1a`). Each is
given as its LaTeX source and as the PDF compiled on the day of the freeze. A frozen text is not edited, so the errata
and open items found since are listed in `ERRATA_NEXT_UNFREEZE.md`, with their places, reasons, and proposed texts. The
folder also holds the design specification v5 behind the design code `reference/CIF_v2.py` and the differences from
the previous freeze (v110). `tests/test_appendix_trace.py` asserts the md5s of the four LaTeX files, and it checks that
`reference/glossary_values.json` holds the preamble's own values, using the preamble's copy in `reference_glossary/`.
`docs/ALIGNMENT.md` lists every statement of the three texts that the code implements, relies on, or must not
contradict.

The code and the documents cite App. A by result and equation, numbered as printed in `APPENDIX_A_revised.pdf` (for
example, Result 8, eq. 47). They cite App. B by claim, premise, or critical question (C4, P-C7.1, CQ5), or by line of
`APPENDIX_B_current.tex` (App. B :4373), the glossary by the name of its subsection, and spec v5 by section (§13).
""")
edit(FR, """| `glossary_official.tex` | `e9dc24fdead748e9df9e9eebd608a62e` (v114: the hash `claude/FREEZE_v114_NOTE.md` names, and the copy `tests/test_appendix_trace.py` asserts; it reads "Depth factor, CFD" for K_depth^CFD, v114's D90.  A `glossary_official.tex` that Review 1 of rc7 reports among the project's uploads, md5 `1c42788c4f7bae9b01d7cf5d54504880` by its account (not among the project's files as attached on 2026-10-21), is by Review 1's diff the copy that reads "measured" there, so it predates v114; Review 1 counts 99 differing lines of wording and citation keys) |
| `Ladder_Array_Design_Spec_v5.md` | `89ecfb3730bbaac632dff60b2483b5bb` |
""", """| `glossary_official.tex` | `e9dc24fdead748e9df9e9eebd608a62e` (v114, which reads "Depth factor, CFD" for K_depth^CFD, v114's D90; an earlier copy reads "measured" there) |
| `Ladder_Array_Design_Spec_v5.md` | `89ecfb3730bbaac632dff60b2483b5bb` |
| `v110_to_v114_APPENDIX_A.diff` | `49b6fa487cc52533e788590122c8365f` |
| `v110_to_v114_APPENDIX_B.diff` | `1c4289ea0d8f17f802fdaacddb90f07c` |
| `v110_to_v114_glossary.diff` | `e610791c85e4b810aba88fa15c972162` |
""")
edit(FR, """The reference kernels the package imports, `opencif/reference/CIF_v2.py` (v2.1.0), are byte-identical to the project's
`CIF_v2.py`: md5 `0e0961a50441afb94671b4470d41bace` (asserted by `tests/test_appendix_trace.py`).
""", """The design code that App. A names, `reference/CIF_v2.py` (v2.1.0), is the project's `CIF_v2.py` byte for byte (md5
`0e0961a50441afb94671b4470d41bace`, asserted by `tests/test_kernels.py` and `tests/test_appendix_trace.py`). The package
does not import it, and the tests compare the package with it.
""")

# ---------------------------------------------------------------------------------------------------------------------
# docs/ALIGNMENT.md
# ---------------------------------------------------------------------------------------------------------------------
AL = 'docs/ALIGNMENT.md'
s = open(os.path.join(TREE, AL), encoding='utf-8').read()
row = next(ln for ln in s.split('\n') if ln.startswith('| AB27 |'))
cells = row.split(' | ')
cells[3] = ('opencif implements Prospectus C7/C8, as a derivative of the design code `CIF_v2.py`; acceptance suite = '
            '`tests/`.  Since 1.0.0 it runs CIF_v2\'s physics auditor, ported (`opencif/checks/auditor.py`), on the '
            'photomask file of every build: P0 (the file is the drawing), R0 (its own ruler on the straight rows), '
            'P1–P6b, P9, R6b and P10 are pass conditions in every mode; the design code\'s M14 is enforced on the drawn '
            'mask (check M14, envelope); P8\'s two criteria (CIF_v2 `_PHYSICS_SRC`: inlet Re ≤ 1.0 at the device\'s '
            'flow, default 5 µL/min; 0.05 Re D_h ≤ Λ/2) are checked on the drawn inlet at the design flow and at 5 µL/min '
            'and reported under their own name (`validation.physics_audit`, check P8), not enforced: at their own flows '
            'every shipped design exceeds the Re limit (3.6–61), at 5 µL/min all pass.  The warrant sentence '
            '(:4375–4377) refers back to the auditor of :4349–4350: it holds of opencif\'s auditor, which enforces M14 and '
            'reports P8.  D95 (up to RC8: the auditor not run on the mask) is closed into D106, revised for 1.0.0, which '
            'names opencif and says which gate its auditor enforces (`docs/DISSERTATION_NOTES.md` § 4); the frozen text '
            'is not edited')
cells[4] = ('trace "P8 (CIF_v2 _PHYSICS_SRC) ...", "P8 inlet Re in the record = CIF_v2 reynolds"; `tests/test_auditor.py` '
            '(the port against CIF_v2\'s own auditor; seeded faults)')
cells[5] = 'ALIGNED (1.0.0; P8 reported, LIMITATION 6); D106 (proposed, revised for 1.0.0) |'
open(os.path.join(TREE, AL), 'w', encoding='utf-8').write(s.replace(row, ' | '.join(cells)))
n += 1
edit(AL, "| AG7i | η (measured, 0.903) | a realised-yield factor | not applied | LIMITATION 5 |",
     "| AG7i | η (measured, 0.903) | a realised-yield factor | not applied; since 1.0.0 ηF is reported beside F on every "
     "record (`validation.design_code_reports.B6_F_real`), as the design code reports it | LIMITATION 5 |")
edit(AL, "Not used: `\\kFval`, `\\ctwobar` (the other pair), `\\etaVal` (LIMITATION 5).",
     "Not used: `\\kFval`, `\\ctwobar` (the other pair).  `\\etaVal` 0.903 = CIF_v2's `Pairing.eta`, which 1.0.0 uses "
     "to report ηF, never applied (LIMITATION 5; the trace test asserts the value against CIF_v2's).")
edit(AL, "9, 11 and 19 fields, so a field added, dropped or left without a disposition fails it; the counts are Pairing 2 gated, 1\npartly gated, 1 reported, 3 consumed, 2 other;",
     "9, 11 and 19 fields, so a field added, dropped or left without a disposition fails it; the counts are Pairing 2 gated, 1\npartly gated, 3 reported, 3 consumed, 0 other;")
edit(AL, "**Reported** (evaluated on every design, never a pass condition): f_soft.",
     "**Reported** (evaluated on every design, never a pass condition): f_soft; and, since 1.0.0, η (ηF, the design "
     "code's F_real) and β_cal (the transfer bias it implies, the design code's mandatory report; on the cubic pair it "
     "does not enter N), each with its value asserted equal to CIF_v2's.")
edit(AL, "its interval, the cubic kernel. **A spec input, replaced, or a named limitation**: the process and geometry inputs; η\n(LIMITATION 5), β_cal (not needed on the cubic pair), μ and ρ_f",
     "its interval, the cubic kernel. **A spec input, replaced, or a named limitation**: the process and geometry inputs; μ and ρ_f")
edit(AL, """M16 against `CIF_v2.march` on matched requests, just below and just above each limit; w_min against `Process.w_min`; M14
against CIF_v2's literal row-1 formulas).  (The trace test's disposition text for L_eq still names only the exact bend
solution; it is a stated, not a tested, disposition, and editing the test would change the tests hash every record
carries, so it waits for the next unfreeze.)""", """M16 against `CIF_v2.march` on matched requests, just below and just above each limit; w_min against `Process.w_min`; M14
against CIF_v2's literal row-1 formulas).  (1.0.0: the trace test's disposition text for L_eq names both of its
replacements, the exact bend solution for arcs and c_corner for square turns; RC8 left it, to keep its tests hash.)""")

# ---------------------------------------------------------------------------------------------------------------------
# docs/DISSERTATION_NOTES.md
# ---------------------------------------------------------------------------------------------------------------------
DN = 'docs/DISSERTATION_NOTES.md'
between(DN, """* **The physics auditor (App. B: "the physics auditor separately checks the drawn device's inlet against its own,""",
        "* **The last rows and the outlet.**",
        """* **The physics auditor (App. B: "the physics auditor separately checks the drawn device's inlet against its own,
  stricter channel limit"; P8).**  Since 1.0.0 opencif runs the design code's physics auditor, ported, on the
  photomask file of every build, read back from the DXF (`validation.physics_audit`, `assessment.status.physics_audit`;
  `docs/PHYSICS_AUDIT.md`).  Say so precisely: the auditor is `CIF_v2.py`'s (P1–P10), ported, with three changes of
  formula a DXF file forces (P1 anchors the cubic pair; P6b uses App. A's march–outlet identity; P9, R6b and P10 are
  measured on the polygon), three checks added (P0, the file is the drawing; R0, an independent ruler on the straight
  rows; R6b, the design code's own acceptance check of the posts, on every post) and one change of role: P8, the
  stricter inlet limit (Re ≤ 1, a settable default; App. B gives no value), is checked on the drawn inlet and recorded
  under its own name, not enforced, because at that limit every lineage device fails at its own flow (Dinh V1's inlet
  Re is 4.65 at 30 µL/min; CIF_v2 runs its lineage demo at 2 µL/min), while every opencif design meets it at 5 µL/min
  (inlet Re 0.38–0.93).  The envelope the thesis claims, the design code's channel gate Re ≤ 10 and the Stokes gates
  (M14), is enforced on the drawn mask (check M14).  The auditor's independence has a stated boundary: its own ruler
  re-measures the straight rows (most of every device's gaps) from the file, and the head units, U-turns and outlet are
  taken from the generator's measurement.  App. B's sentence that the software is "warranted by their acceptance suite
  and the physics auditor" then holds of opencif's auditor, with the difference that it reports rather than enforces
  the stricter limit; the erratum states it (D106, revised for 1.0.0, into which D95 is closed;
  `foundation/ERRATA_NEXT_UNFREEZE.md` in the dissertation's copy):
  > A:184–187, after "the design code `CIF_v2.py` (Results 10–14)": "whose kernels the targets-to-photomask generator
  > opencif 1.0.0 (code stamp `{{STAMP}}`) reproduces bit for bit"; B:4375–4377: "... warranted by their acceptance
  > suite and the physics auditor (the design code's, which opencif 1.0.0 runs on every photomask file: it enforces the
  > design code's channel and Stokes gates on the drawn mask and reports its own stricter inlet limit) rather than by a
  > defeasible inference."
  Do not call `CIF_v2.py` "the reference implementation": in the frozen texts that phrase names the preprint's
  notebook; App. A's name for it is "the design code".  Cite the port's validation (`tests/test_auditor.py`: the same
  numbers as CIF_v2's own auditor on its demonstration devices; a seeded fault for every check) and the release's audit
  of every shipped photomask (`docs/VALIDATION.md`).
""")

edit(DN, """     taper-direction sign test (Fig. 4).

## 4. Claims to make carefully""", """     taper-direction sign test (Fig. 4).

8. **1.0.0: the physics audit of every photomask, with its own ruler, and the design code's reports.**  The
   design code's physics auditor (`CIF_v2.py`'s P1–P10), ported and run on the photomask file of every build, read back
   from the DXF (`opencif/checks/auditor.py`; `docs/PHYSICS_AUDIT.md`), with three checks added: P0 (the file is the
   drawing), R0 (an independent ruler, `opencif/checks/ruler.py`, re-measures every ordinary straight row from the file
   and holds the generator's measurement to it within 10⁻³) and R6b (the design code's acceptance check of the posts, on
   every post); the design code's gate M14 enforced on the drawn mask; and the design code's mandatory reports and
   acceptance validators (pressure drop, ηF, transfer bias, sensitivity map, porosity, operating regime, steric screen,
   fabrication tolerance) recorded on every build by its own formulas (`opencif/checks/design_code.py`).  Evidence: the
   port against CIF_v2's own auditor and seeded faults (`tests/test_auditor.py`), each report against CIF_v2's own
   function (`tests/test_design_code.py`), the kernels bit for bit (`tests/test_kernels.py`), and the audit of every
   shipped photomask, in the build and again from scratch (`docs/VALIDATION.md`).  `docs/CIF_V2_REQUIREMENTS.md`
   lists every requirement of `CIF_v2.py` and its status.

## 4. Claims to make carefully""")

# ---------------------------------------------------------------------------------------------------------------------
# RC8's module names, in the documents that describe the current code
# ---------------------------------------------------------------------------------------------------------------------
KERNEL_MAP = {'g_cubic': 'g_cubic', 'r': 'r', 'x_of_f': 'x_of_f'}


def rename(rel):
    global n
    p = os.path.join(TREE, rel)
    s = open(p, encoding='utf-8').read()
    s0 = s
    # "`cifsched.f_of` → `CIF_v2.g_cubic`": the schedule now calls the restated kernels, equal to CIF_v2's
    s = re.sub(r'`cifsched\.(\w+)` → `CIF_v2\.(\w+)(\([^`]*\))?`',
               lambda m: '`schedule.%s` → `kernels.%s%s` (= `CIF_v2.%s`, bit for bit)' % (m.group(1), m.group(2),
                                                                                       m.group(3) or '', m.group(2)), s)
    s = re.sub(r'`cifsched\.(\w+)`', r'`schedule.\1`', s)
    s = s.replace('`cifsched.py`', '`opencif/physics/schedule.py`').replace('`measure2.py`', '`opencif/checks/measure.py`')
    s = re.sub(r'`measure2\.(\w+)`', r'`measure.\1`', s)
    s = re.sub(r'`cifsched`', '`schedule`', s)
    if s != s0:
        open(p, 'w', encoding='utf-8').write(s)
        n += 1


for rel in ('docs/ALIGNMENT.md', 'docs/DISSERTATION_NOTES.md', 'docs/LIMITATIONS.md', 'docs/PRESET_PROVENANCE.md'):
    rename(rel)
edit(DN, "`reference/CIF_v2.py` (byte-frozen kernels), `schedule.f_of`, `x_of`, `KAPPA_G`",
     "`opencif/physics/kernels.py` (the design code's kernels, restated and tested bit for bit against the frozen "
     "`reference/CIF_v2.py`), `schedule.f_of`, `x_of`, `KAPPA_G`")

# ---------------------------------------------------------------------------------------------------------------------
# the README: RC8's becomes the reference, a new one opens the repository; the new documents
# ---------------------------------------------------------------------------------------------------------------------
shutil.move(os.path.join(TREE, 'README.md'), os.path.join(TREE, 'docs', 'REFERENCE.md'))
for f in ('README.md',):
    shutil.copy(os.path.join(HERE, 'docs', f), os.path.join(TREE, f))
for f in ('ARCHITECTURE.md', 'PHYSICS_AUDIT.md', 'VALIDATION.md', 'CIF_V2_REQUIREMENTS.md', 'EVALUATION.md'):
    shutil.copy(os.path.join(HERE, 'docs', f), os.path.join(TREE, 'docs', f))
RF = 'docs/REFERENCE.md'
edit(RF, "# opencif\n\nTarget-first generator", "# opencif: the reference\n\nThe full reference of the generator (RC8's README, kept current for 1.0.0; `README.md` is the starting\npage).  Target-first generator")
edit(RF, "(`cifsched.bend_ratio`)", "(`schedule.bend_ratio`)")
edit(RF, "(`cifsched.outlet_exp_c`;", "(`schedule.outlet_exp_c`;")
edit(RF, """  - `assessment.status` gives the verdict's parts as separate fields (hard checks, envelope, chip area, aspect ratio,
    G9 length budget, the drawn-mask assessment that was run, CIF_v2's physics auditor (not run on the mask; P8's
    criteria are diagnostics), that nothing was fabricated or measured, and, for a reproduction, what the label
    covers).""", """  - `assessment.status` gives the verdict's parts as separate fields (hard checks, envelope, chip area, aspect ratio,
    G9 length budget, the drawn-mask assessment that was run, the physics audit of the photomask file
    (`status.physics_audit`: PASS, PASS_OUTSIDE_ENVELOPE, FAIL or INCOMPLETE), that nothing was fabricated or
    measured, and, for a reproduction, what the label covers).
* **The physics audit (1.0.0), on the photomask file read back:** the design code's physics auditor (CIF_v2's P1–P10),
  ported, with P0 (the file is the drawing), R0 (the auditor's own ruler re-measures every ordinary straight row from
  the file; the generator's measurement must agree within 10⁻³), R6b (every straight-row post is the process's post)
  and the design code's gate M14 on the drawn mask.  Its hard checks are pass conditions in every mode (exit 2), M14
  is part of the validity envelope (exit 3), P7 and P8 (CIF_v2's stricter inlet limit) are reported
  (`validation.physics_audit`; `docs/PHYSICS_AUDIT.md`).  `python -m opencif audit RUN_DIR` repeats it from scratch.
* **The design code's reports (1.0.0):** the mandatory reports of CIF_v2's `design()` and the operating validators of
  its acceptance run, by its own formulas, recorded, never a pass condition (`validation.design_code_reports`).""")
edit(RF, """* **Parts taken from the design, not measured:**
  - the post length across the flow (gap length l_g), and so the side channels' inner edge;
  - the radii of the head and outlet arcs (their widths are measured);
  - the head rows' centre-path lengths;
  - any off-chip tubing.""", """* **Parts taken from the design, not measured:**
  - the radii of the head, U-turn and outlet arcs (their widths are measured);
  - the head rows' centre-path lengths;
  - any off-chip tubing.

  The straight rows' lengths and the post length across the flow (l_g) are measured on the file by the physics
  audit's own ruler (R0) and its post check (R6b) since 1.0.0 (LIMITATIONS 12).""")
edit(RF, """python3 tests/test_units.py                                                  # unit checks (a few minutes)""",
     """python3 -m opencif audit runs/<name>                                        # the physics audit again, from the photomask file (1.0.0)
python3 tests/test_kernels.py                                                # (1.0.0) the kernels = CIF_v2.py's, bit for bit
python3 tests/test_auditor.py                                                # (1.0.0) the auditor against CIF_v2's own; seeded faults
python3 tests/test_design_code.py                                            # (1.0.0) the design code's reports against its own functions
python3 tests/test_teaching.py                                               # (1.0.0) the teaching module = the generator's plan
python3 tests/test_units.py                                                  # unit checks (a few minutes)""")
edit(RF, """#   (needs the element's runs/<name>/<name>.pkl, which the evaluation archive does not ship; it writes runs/<name>_x3/)""",
     """#   (needs the element's runs/<name>/<name>.pkl from a build on this machine; it writes runs/<name>_x3/)""")
edit(RF, """| 2 | FAIL; the photomask is kept for diagnosis as `<name>_FAILED.dxf`, never under the deliverable name `<name>.dxf` (rc7) |""",
     """| 2 | FAIL (a hard check, the physics audit's included); the photomask is kept for diagnosis as `<name>_FAILED.dxf`, never under the deliverable name `<name>.dxf` (rc7) |""")
os.remove(os.path.join(TREE, 'EVALUATION.md'))                   # one copy: docs/EVALUATION.md
# the changelog: 1.0.0 above RC8
edit('docs/CHANGELOG.md', '## 1.0.0rc8 (freeze 2026-10-21)',
     open(os.path.join(HERE, 'docs', 'CHANGELOG_1.0.0.md'), encoding='utf-8').read() + '## 1.0.0rc8 (freeze 2026-10-21)')
# no em dashes in the release's documents: the earlier releases' entries keep their words otherwise, an empty cell of
# their tables reads "none", and their one unnumbered item reads "unnumbered".
edit('docs/CHANGELOG.md', '\n| \u2014 | `validation.gates.plateau_u_min_all_planned_rows`',
     '\n| unnumbered | `validation.gates.plateau_u_min_all_planned_rows`')
edit('docs/CHANGELOG.md', '| \u2014 |', '| none |', count=11)
for old, new in (
        ("the outlet stretch-end term's domain \u2014 every record states it\n",
         "the outlet stretch-end term's domain, which every record states\n"),
        ("`gifford_plt()` removed \u2014\n  they carried", "`gifford_plt()` removed,\n  because they carried"),
        ('`tests/test_negative.py`: 12b \u2014 `reproduction.L_out`', '`tests/test_negative.py`: 12b: `reproduction.L_out`'),
        ('60 µm deep, the narrow head) \u2014 an explicit', '60 µm deep, the narrow head): an explicit'),
        ('67); 13 \u2014 the stretch-end domain', '67); 13: the stretch-end domain'),
        ("the shipped masks' panels \u2014 Dinh p1–p4, Gifford Fig. 4 and\n  Fig. 6, Xia c–g and Fig. 3a \u2014 each",
         "the shipped masks' panels (Dinh p1–p4, Gifford Fig. 4 and\n  Fig. 6, Xia c–g and Fig. 3a), each")):
    edit('docs/CHANGELOG.md', old, new)
edit('requirements-tested.txt', '# the exact versions the evidence in runs/ was produced with (Python 3.11.15); requirements.txt gives the minimums',
     '# the exact versions every record of opencif 1.0.0 was produced with (Python 3.11.17); requirements.txt gives the\n'
     '# minimums.  Under Python 3.13.16 with these versions, an independent reviewer\'s run of RC8\'s tutorial gave the same\n'
     '# printed numbers and the same photomask bytes.')
edit('requirements.txt', '# opencif 1.0.0rc1 -- tested with Python 3.11.15 and the versions below',
     '# opencif 1.0.0 -- tested with Python 3.11.17 and the versions below')
# the repository's own files: citation, packaging, CI, git attributes
for rel in ('LICENSE', 'CITATION.cff', 'pyproject.toml', '.gitignore', '.gitattributes', '.github/workflows/tests.yml'):
    dst = os.path.join(TREE, *rel.split('/'))
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy(os.path.join(HERE, 'meta', *rel.split('/')), dst)
# RC8's discrepancy register, which the tutorial and LIMITATIONS cite, is kept with the history; the sketch of Figure 1
# (README and tutorial) is drawn by a tool of its own
shutil.copy(os.path.join(HERE, 'meta', 'provenance', 'RC8_DISCREPANCIES.md'),
            os.path.join(TREE, 'provenance', 'RC8_DISCREPANCIES.md'))
shutil.copy(os.path.join(HERE, 'repo_tools', 'tutorial_figures.py'), os.path.join(TREE, 'tools', 'tutorial_figures.py'))
print('patch_docs: %d edits' % n)
