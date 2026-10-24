"""Early loss of beads of diameter d along a built CIF device (my own re-implementation of the INV-2 streamline idea):
beads start spread evenly in flux coordinate psi over [psi_min, 1 - psi_min], psi_min = the flux fraction at which the
record's realised cutoff equals d (steric exclusion); row i extracts psi < fA(i) at side A and psi > 1 - fB(i) at side B;
the rest is remapped psi -> (psi - fA)/(1 - fA - fB) and held at least psi_min from either wall.  Prints the share of beads extracted in the first rows and the
expected count in leg 4's outer band (INV-2: the outer 0.2 holds the earliest 13.5 % of the inlet flow and would hold
>= 98 beads at the inlet concentration; 2 are seen).
usage: python early_loss.py TREE PKL d_um [rows]"""
import pickle, sys, math
import numpy as np
T = sys.argv[1]
sys.path[:0] = [T, T + '/opencif', T + '/opencif/engine']
a = pickle.load(open(sys.argv[2], 'rb'))['arrays']
d = float(sys.argv[3]); M = int(sys.argv[4]) if len(sys.argv) > 4 else 300
fA, fB, dA, dB = (np.asarray(a[k], float) for k in ('fA', 'fB', 'dA', 'dB'))
D = np.concatenate([dA, dB]); F = np.concatenate([fA, fB]); ok = np.isfinite(D) & (F > 0)
D, F = D[ok], F[ok]; o = np.argsort(D); D, F = D[o], np.maximum.accumulate(F[o])
psi_min = float(np.exp(np.interp(math.log(d), np.log(D), np.log(F))))
psi = np.linspace(psi_min, 1 - psi_min, 200001)[1:-1]
alive = np.ones(psi.size, bool); lost = 0
rows_above = [i + 1 for i in range(len(dA)) if max(dA[i], dB[i]) > d]
for i in range(min(M, len(fA))):
    out = alive & ((psi < fA[i]) | (psi > 1 - fB[i]))
    lost += out.sum(); alive &= ~out
    psi = np.where(alive, np.clip((psi - fA[i]) / (1 - fA[i] - fB[i]), psi_min, 1 - psi_min), psi)   # steric: never nearer a wall than psi_min
share = lost / psi.size
print('d %.1f um: rows with a realised cutoff above d: %s; beads lost in the first %d rows: %.3f %%; '
      'expected in the outer band ~%.1f (seen 2)' % (d, rows_above[:12] if len(rows_above) <= 12 else '%d rows (%d..%d)' % (len(rows_above), rows_above[0], rows_above[-1]),
                                                   M, 100 * share, 98 * share / 0.135))
print('row 1-5 realised cutoffs (um):', np.round(np.maximum(dA[:5], dB[:5]), 2), ' declared rows above D:', int(np.sum(np.maximum(dA, dB) > 1.01 * np.nanmedian(dA[200:]))))
