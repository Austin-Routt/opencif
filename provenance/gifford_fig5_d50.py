"""This project's reading of Gifford 2014 Fig. 5 (the held-out bead devices A-D): the fraction of beads of each size each
long device retains, at its design cutoff D (the preset's targets.cutoff_um) and as a 50 %-retention size d50.

Input: the enrichment E of each bead size (4, 6, 8, 10 um) in the retentate, as printed in the Fig. 5 caption (p. 4501).
Two values are assumed, not printed: 5C "could no longer enrich 4 um beads" is read as E(4) = 1, and 5D prints only
"enrichment of 12.0x for the largest 10 um beads", so its 4, 6 and 8 um beads are taken as not enriched (E = 1).  Any
E(8) <= 1 gives R(8) <= 1 - F and moves 5D's d50 by under 0.02 um; neither assumption moves R(D) of 5A-5C.
A bead that is fully retained is enriched by 1/(1 - F), so the retained fraction is R = E (1 - F).  Three values of the
ceiling 1/(1 - F) are carried, because the result depends on it:
  duct law  1 - F from the outlet widths (w_c 100 um, w_s 450 um, depth 125 um) through the duct law (App. A R14): 23.5x;
  caption   "equivalent to a theoretical 24.4x enrichment" (Fig. 5 caption);
  text      "theoretically ~23 : 1" side : centre flow (p. 4500), i.e. 24x.
The duct-law ceiling makes R(10 um) exceed 1 for 5A and 5B (24.1x and 24.2x measured): the data themselves put the
ceiling at 24.2x or more, which favours the caption's 24.4x.

What this is and is not.  kappa_g was fitted on Gifford's Fig. 3 thresholds, "the filtration fraction per gap below which
the beads of a specific diameter are consistently maintained in the central flow channel" (p. 4500), measured on
"short-length (several centimeter long) test devices" (p. 4499): an onset-of-loss criterion on a short array.  Fig. 5's
devices are long (1575-7161 gaps here) and report cumulative retention.  So R(D) and d50 measure how much a long device
loses at and near D; they are not a test of kappa_g at its +-3.4 % level (the adjudication review's E3: a d50-within-CI
pass test is ill-posed under an onset criterion).  Gifford: "there are invariably places in a device where some fraction
of these beads will be lost ... often owing to imperfections ... which can serve to elevate the effective value of fgap"
(p. 4500), and "when a value of fgap is chosen that generates an array with a cutoff size several microns below a given
particle size, ... the vast majority of those particles are indeed retained" (p. 4500).
New evidence of this project: not part of App. A, App. B or the glossary.
Run from the repository root: python3 provenance/gifford_fig5_d50.py"""
import math, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import yaml
from opencif import build  # noqa: F401
from opencif.physics import schedule as S

E = {'A': {4: 8.6, 6: 17.7, 8: 22.2, 10: 24.1}, 'B': {4: 2.3, 6: 16.9, 8: 22.6, 10: 24.2},
     'C': {4: 1.0, 6: 2.3, 8: 19.3, 10: 22.7}, 'D': {4: 1.0, 6: 1.0, 8: 1.0, 10: 12.0}}
# 5C: "could no longer enrich 4 um beads" -> E(4) = 1; 5D: only the 10 um beads enriched -> E(4, 6, 8) = 1 (assumed)
CEIL = {'duct law': 1 - 2 * S.rr(100, 125) / (S.rr(450, 125) + 2 * S.rr(100, 125)),
        'caption 24.4x': 1 / 24.4, 'text ~23:1': 1 / 24.0}
HOW = ('lin', 'logd', 'logit')


def _lg(r):
    r = min(max(r, 1e-6), 1 - 1e-6)
    return math.log(r / (1 - r))


def at(R, x, how):
    """R interpolated at bead size x"""
    ds = sorted(R)
    for x0, x1 in zip(ds[:-1], ds[1:]):
        if x0 <= x <= x1:
            if how == 'lin':
                return R[x0] + (x - x0) / (x1 - x0) * (R[x1] - R[x0])
            if how == 'logd':
                return R[x0] + (math.log(x) - math.log(x0)) / (math.log(x1) - math.log(x0)) * (R[x1] - R[x0])
            z = _lg(R[x0]) + (x - x0) / (x1 - x0) * (_lg(R[x1]) - _lg(R[x0]))
            return 1 / (1 + math.exp(-z))
    return None


def d50(R, how):
    ds = sorted(R)
    for x0, x1 in zip(ds[:-1], ds[1:]):
        if R[x0] < 0.5 <= R[x1]:
            if how == 'lin':
                return x0 + (x1 - x0) * (0.5 - R[x0]) / (R[x1] - R[x0])
            if how == 'logd':
                t = (0.5 - R[x0]) / (R[x1] - R[x0])
                return math.exp(math.log(x0) + t * (math.log(x1) - math.log(x0)))
            return x0 + (0 - _lg(R[x0])) / (_lg(R[x1]) - _lg(R[x0])) * (x1 - x0)
    return None


Dd = {dev: yaml.safe_load(open(os.path.join(ROOT, 'presets', 'gifford2014_fig5%s_as_fabricated_12p5um.yaml' % dev)))
      ['targets']['cutoff_um'] for dev in E}
summary = {dev: dict(RD=[], q=[], above=False) for dev in E}
for cn, c in CEIL.items():
    print('ceiling %-13s 1 - F = %.4f (full retention = %.1fx)' % (cn, c, 1 / c))
    for dev, e in E.items():
        R = {k: float(v * c) for k, v in e.items()}
        rD = {h: at(R, Dd[dev], h) for h in HOW}
        dd = {h: d50(R, h) for h in HOW}
        summary[dev]['RD'] += list(rD.values())
        summary[dev]['q'] += [v / Dd[dev] for v in dd.values() if v]
        summary[dev]['above'] |= any(v is None for v in dd.values())
        print('  Fig. 5%s  D %.3f um  R %s%s  R(D) %s  d50/D %s' % (
            dev, Dd[dev], {k: round(v, 3) for k, v in R.items()}, '  (R > 1: ceiling too low)' if max(R.values()) > 1 else '',
            '/'.join('%.3f' % rD[h] for h in HOW), '/'.join(('%.3f' % (dd[h] / Dd[dev])) if dd[h] else '>10um' for h in HOW)))
print('summary over the three ceilings and three interpolants (lin / log d / logit):')
for dev, v in summary.items():
    print('  Fig. 5%s  D %.3f um: R(D) %.2f-%.2f; d50/D %.3f-%.3f%s' % (
        dev, Dd[dev], min(v['RD']), max(v['RD']), min(v['q']), max(v['q']),
        ' (above 10 um, i.e. d50/D > 1.108, under the caption ceiling)' if v['above'] else ''))
