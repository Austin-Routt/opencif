"""rc6: the fixes from the independent review of rc6's code delta (2026-10-17, provenance/RC6_REVIEWS_VERIFICATION_
2026-10-17.md section 14.2).  Each edit by exact string replacement (it stops if a target is not found exactly once).

  R5   engine/devspecs.py: xia_fig's docstring is formatted with the rc6 scale only when it exists (under python -OO
       docstrings are stripped and the rc6 line raised TypeError on import; rc5 imported fine).
  R7   build.py: validation.features_found also counts the missing features by kind (by_kind), so that a drawing that
       loses more than the 20 listed entries still says which kinds it lost (the list keeps the first 20).
  R11  pipeline.py: the strict refusal of an outlet stretch-end beyond the solved tables names the outlet style, as every
       other outlet refusal does (tests/test_lattice.py requires it of every exit-5 outlet refusal; latent: no lattice
       case reaches it).
  R12  engine/cifsched.py: the comments of the ring defaults say they are rc5's Fig. 3a reading and unused by the
       shipped specs (spec.to_engine always sets them); the values are unchanged.

Usage: python3 apply_rc6_review_fixes.py OPENCIF_ROOT"""
import os
import sys

ROOT = sys.argv[1]


def edit(rel, pairs):
    p = os.path.join(ROOT, rel)
    s = open(p).read()
    for old, new in pairs:
        n = s.count(old)
        if n != 1:
            raise SystemExit('%s: target found %d times: %r' % (rel, n, old[:80]))
        s = s.replace(old, new)
    open(p, 'w').write(s)
    print('edited', rel, len(pairs))


edit('opencif/engine/devspecs.py', [
    ("xia_fig.__doc__ = xia_fig.__doc__ % R6['xia']['scale_um_per_px']     # rc6: the scale, from rc6_values.json\n",
     "if xia_fig.__doc__:                    # rc6: the scale, from rc6_values.json (no docstring under python -OO)\n"
     "    xia_fig.__doc__ = xia_fig.__doc__ % R6['xia']['scale_um_per_px']\n"),
])
edit('opencif/build.py', [
    ("""    V['features_found'] = dict(ok=not m['not_found'], not_found=m['not_found'][:20], n_not_found=len(m['not_found']))""",
     """    V['features_found'] = dict(ok=not m['not_found'], not_found=m['not_found'][:20], n_not_found=len(m['not_found']),
                               by_kind={kd: sum(1 for q in m['not_found'] if q['kind'] == kd)
                                        for kd in sorted({q['kind'] for q in m['not_found']})})"""),
])
edit('opencif/pipeline.py', [
    ("""            opts['outlet_refusal'] = ('%s: build it in evaluate or reproduce mode, where the end values are held and '
                                      'reported, or change the outlet' % e)""",
     """            opts['outlet_refusal'] = ('the %s outlet: %s: build it in evaluate or reproduce mode, where the end values '
                                      'are held and reported, or change the outlet' % (opts.get('outlet'), e))"""),
    ("""        opts['outlet_refusal'] = ('%s: build it in evaluate or reproduce mode, where the end values are held and '
                                  'reported, or change the outlet' % e)""",
     """        opts['outlet_refusal'] = ('the %s outlet: %s: build it in evaluate or reproduce mode, where the end values '
                                  'are held and reported, or change the outlet' % (opts.get('outlet'), e))"""),
])
edit('opencif/engine/cifsched.py', [
    ("""    pc_inlet_offset: float = 0.0     # paperclip: the inlet port lies this far (um) downstream of the loop-2 centres
                                     # (Xia Fig. 3a: 2.54 mm)""",
     """    pc_inlet_offset: float = 0.0     # paperclip: the inlet port lies this far (um) downstream of the loop-2 centres
                                     # (Xia Fig. 3a: 2.55 mm at rc6's reading; the preset sets it)"""),
    ("""    pc_out_ext: Optional[float] = None  # paperclip, ring outlet: the ring's top lies this far (um) beyond the centres of
                                     # the loops at the last leg's end (Xia Fig. 3a: 3.43 mm); None = clearance minimum""",
     """    pc_out_ext: Optional[float] = None  # paperclip, ring outlet: the ring's top lies this far (um) beyond the centres of
                                     # the loops at the last leg's end (Xia Fig. 3a: 3.44 mm at rc6's reading; the preset
                                     # sets it); None = clearance minimum"""),
    ("""    ring_R: float = 4179.0           # 'ring': centreline radius of the circle (um)""",
     """    ring_R: float = 4179.0           # 'ring': centreline radius of the circle (um).  These three defaults are rc5's
                                     # Fig. 3a reading; no shipped spec uses them (spec.to_engine sets every ring value;
                                     # rc6's Xia values are in data/r6/rc6_values.json)"""),
])
