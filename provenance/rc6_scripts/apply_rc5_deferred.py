"""rc6: apply the edits the rc5 verification deferred to rc6 (provenance/RC5_REVIEWS_VERIFICATION_2026-10-12.md
sec. 8.2), in place, each by exact string replacement (it stops if a target is not found exactly once):

 (1) the stretch-end tables' domain: a status in every record (design.outlet.stretch_end.table_domain), a refusal in
     strict mode (exit 5, 'REFUSED (outlet)'), and a clause in the verdict in evaluate and reproduce modes.  The edits
     are the report's, as tried on a scratch copy of rc5 (provenance/rc5_scripts/rc6_patch1_trial/apply_patch1.py),
     plus the verdict clause the trial left out;
 (2) reproduction.L_out with outlet.flare_deg refused (a SpecError naming both);
 (3) the fail-safe refusal's wording.

(4), the lead arc's interaction in c for teardrop leads, is not applied (optional in the report; it moves the nine
teardrop designs' last rows by a fraction of their model band, which already contains it).  (5) is a test-only
option, not applied.

Usage: python3 apply_rc5_deferred.py OPENCIF_ROOT"""
import os
import sys

R = sys.argv[1]


def edit(rel, old, new, count=1):
    p = os.path.join(R, rel)
    s = open(p).read()
    assert s.count(old) == count, (rel, old[:70], s.count(old))
    open(p, 'w').write(s.replace(old, new))


# ---- (1) cifsched.py: the domain status, after outlet_exp_len
edit('opencif/engine/cifsched.py', '''def outlet_exp_len(spec, w_a, w_b):
    """the same excess as a length (um) of the narrower channel"""
    return outlet_exp_c(spec, w_a, w_b) * min(w_a, w_b)
''', '''def outlet_exp_len(spec, w_a, w_b):
    """the same excess as a length (um) of the narrower channel"""
    return outlet_exp_c(spec, w_a, w_b) * min(w_a, w_b)


def outlet_exp_domain(spec, w_a, w_b):
    """None where the stretch-end term is inside the solved 2-D tables; else what was held at a table's end (rc6;
    a review of rc5: np.interp clamps silently).  Steps: deep tabulated to R 8 (the shallow closed form has no end);
    flares: R 4 and 45 deg (flare_deg is limited to [1, 45] by the spec)."""
    wn, ww = min(w_a, w_b), max(w_a, w_b)
    if getattr(spec, 'junction_band', None) == 'none' or ww <= wn * (1.0 + 1e-6):
        return None
    R = ww / wn
    if getattr(spec, 'flare_deg', None) is not None and w_b > w_a:
        return ('flare width ratio %.2f beyond the tables\\' %.1f: deep and shallow held there' % (R, FLARE_RATIOS[-1])
                if R > FLARE_RATIOS[-1] + 1e-9 else None)
    tab = STEP_ONE_DEEP if w_b > w_a else STEP_SYM_DEEP
    return ('step width ratio %.2f beyond the deep table\\'s %.1f: deep held there (shallow: closed form)' % (R, tab[-1][0])
            if R > tab[-1][0] + 1e-9 else None)
''')

# spec.py: the strict-mode switch reaches the engine Spec
edit('opencif/spec.py', '''    spec = S.Spec(**kw)
''', '''    spec = S.Spec(**kw)
    # rc6 (rc5 review P, deferred): in strict mode an outlet stretch-end beyond the solved 2-D tables is refused
    spec.exp_domain_refuse = (g(d, 'gates', 'mode', 'strict') == 'strict')
''')

# layout.py: the exception and the checks
edit('opencif/engine/layout.py', '''class OutletInfeasible(ValueError):''', '''class OutletOutsideTables(Exception):
    """rc6: the outlet stretch's end lies outside the solved 2-D tables in strict mode.  A plain Exception, so the leg
    fitter treats it neither as an OutletInfeasible ('try the next leg count') nor as a ValueError ('n legs not
    drawable'): strict mode refuses at the first leg count whose settled lead lies beyond a table.  That is the
    conservative choice; another count could move the lead width, through the stretch length."""


def check_exp_domain(spec, w_last, w_lead):
    dom = S.outlet_exp_domain(spec, w_last, w_lead)
    if dom and getattr(spec, 'exp_domain_refuse', False):
        raise OutletOutsideTables('the outlet stretch\\'s end lies outside the solved 2-D tables (%s)' % dom)
    return dom


def check_plan_domain(sc2):
    dom = sc2.info.get('exp_domain')
    if dom and getattr(sc2.spec, 'exp_domain_refuse', False):
        raise OutletOutsideTables('the outlet stretch\\'s end lies outside the solved 2-D tables (%s)' % dom)


class OutletInfeasible(ValueError):''')
# the two plan-stage placements, after the final outlet_Rc_lead call, on the widths that call used
edit('opencif/engine/layout.py', '''    Rc_lead, Lc, fr = outlet_Rc_lead(S.reclose(sc, turns, R_turn, spec.R_ext_c))
    sc2 = S.reclose(sc, turns, R_turn, Rc_lead + spec.R_ext_c)''', '''    sc_c = S.reclose(sc, turns, R_turn, spec.R_ext_c)
    Rc_lead, Lc, fr = outlet_Rc_lead(sc_c)
    dom = S.outlet_exp_domain(sc_c.spec, sc_c.ws[sc_c.N - 2], float(fr['w'] if 'w' in fr else fr['wl']))
    sc2 = S.reclose(sc, turns, R_turn, Rc_lead + spec.R_ext_c)
    sc2.info['exp_domain'] = dom                              # rc6: checked in plan_layout on the accepted plan''')
edit('opencif/engine/layout.py', '''    Rc_lead, Lc, fr = outlet_Rc_lead(S.reclose(sc, turns, Rt, sc.spec.R_ext_c))
    sc2 = S.reclose(sc, turns, Rt, Rc_lead + sc.spec.R_ext_c)''', '''    sc_c = S.reclose(sc, turns, Rt, sc.spec.R_ext_c)
    Rc_lead, Lc, fr = outlet_Rc_lead(sc_c)
    dom = S.outlet_exp_domain(sc_c.spec, sc_c.ws[sc_c.N - 2], float(fr['w'] if 'w' in fr else fr['wl']))
    sc2 = S.reclose(sc, turns, Rt, Rc_lead + sc.spec.R_ext_c)
    sc2.info['exp_domain'] = dom                              # rc6: checked in plan_layout on the accepted plan''')
# plan_layout checks only the plan its loop-radius search accepts
edit('opencif/engine/layout.py', '''        if ok:
            return sc2, legs
        R = (sc2.info['leg_d'] if kind == 'paperclip' else sc2.info['R_turn']) * 1.1''', '''        if ok:
            check_plan_domain(sc2)                            # rc6: only the plan the radius search accepts
            return sc2, legs
        R = (sc2.info['leg_d'] if kind == 'paperclip' else sc2.info['R_turn']) * 1.1''')
# the three drawing-time planners, before their exp_info calls
edit('opencif/engine/layout.py', '''    L_lead = side_lead_path(fr, x_e, Lr).length() - 2.0          # to the port edge
    info = dict(kind='separate',''', '''    L_lead = side_lead_path(fr, x_e, Lr).length() - 2.0          # to the port edge
    check_exp_domain(spec, w_last, fr['wl'])                     # rc6: the stretch-end tables' domain
    info = dict(kind='separate',''')
edit('opencif/engine/layout.py', '''    exp_term = exp_info(spec, w_last, g['w'], h)''', '''    check_exp_domain(spec, w_last, g['w'])                       # rc6: the stretch-end tables' domain
    exp_term = exp_info(spec, w_last, g['w'], h)''')
edit('opencif/engine/layout.py', '''    L_lead = ring2_path(g, x_e).length() - 2.0          # minus the 1 um overlaps at both ends
    info = dict(kind='combined', style='ring',''', '''    L_lead = ring2_path(g, x_e).length() - 2.0          # minus the 1 um overlaps at both ends
    check_exp_domain(spec, w_last, w)                            # rc6: the stretch-end tables' domain
    info = dict(kind='combined', style='ring',''')
# exp_info records the domain
edit('opencif/engine/layout.py', '''             excess_len_um=float(c * wn), R=float(S.rr(wn, h) * c * wn * S.UM))''', '''             excess_len_um=float(c * wn), R=float(S.rr(wn, h) * c * wn * S.UM),
             table_domain=S.outlet_exp_domain(spec, w_last, w_lead) or 'inside the solved tables')''')

# pipeline.py: the pre-check and the build handler refuse with their own text; the verdict clause
edit('opencif/pipeline.py', '''        except L.OutletInfeasible as e:
            errs.append('%d leg(s): %s' % (n, e))''', '''        except L.OutletOutsideTables as e:                # rc6: the stretch-end tables' domain (strict mode)
            opts['outlet_refusal'] = ('%s: build it in evaluate or reproduce mode, where the end values are held and '
                                      'reported, or change the outlet' % e)
            say(opts['outlet_refusal'])
            return
        except L.OutletInfeasible as e:
            errs.append('%d leg(s): %s' % (n, e))''')
edit('opencif/pipeline.py', '''    except L.OutletInfeasible as e:                       # found only while fitting or drawing: the same refusal''',
     '''    except L.OutletOutsideTables as e:                    # rc6: the stretch-end tables' domain (strict mode)
        opts['outlet_refusal'] = ('%s: build it in evaluate or reproduce mode, where the end values are held and '
                                  'reported, or change the outlet' % e)
        return _write_refusal(d, name, outdir, spec, sc, g0, notes, opts, env_bad, say)
    except L.OutletInfeasible as e:                       # found only while fitting or drawing: the same refusal''')
edit('opencif/pipeline.py', '''    if model_txt and code in (0, 3):
        verdict += ' [1-D model: %s]' % '; '.join(model_txt)''', '''    # rc6: an outlet stretch-end beyond the solved 2-D tables (refused in strict mode; here held and reported)
    td_ = (((rep.get('design') or {}).get('outlet') or {}).get('stretch_end') or {}).get('table_domain')
    if td_ and td_ != 'inside the solved tables':
        model_txt.append('stretch-end term outside the solved tables: %s' % td_)
    if model_txt and code in (0, 3):
        verdict += ' [1-D model: %s]' % '; '.join(model_txt)''')

# ---- (2) reproduction.L_out with a flare
edit('opencif/spec.py', '''    bad = set(repro) - REPRO_FIELDS
    if bad:
        raise SpecError('reproduction: unknown engine fields %s' % sorted(bad))
''', '''    bad = set(repro) - REPRO_FIELDS
    if bad:
        raise SpecError('reproduction: unknown engine fields %s' % sorted(bad))
    if 'L_out' in repro and kw.get('flare_deg') is not None:    # rc6 (rc5 verification, deferred)
        raise SpecError('reproduction.L_out with outlet.flare_deg: the flare sets the outlet stretch '
                        '(outlet.flare_start_um straight, then the flare and its end arc); give flare_start_um instead')
''')

# ---- (3) the fail-safe refusal's wording
edit('opencif/spec.py', """            raise SpecError('outlet.flare_deg with physics.fail_safe_margin is not supported: the fail-safe re-solve '
                            'moves the lead width the flare\\'s length was planned on')""",
     """            raise SpecError('outlet.flare_deg with physics.fail_safe_margin is not supported: the fail-safe re-solve '
                            'rebuilds the march, which moves the last side channel\\'s width the flare was planned on')""")
print('applied to', R)
