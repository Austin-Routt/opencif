"""package.py: package the code base and the run records for delivery, with every zip under LIMIT_MB.

Usage: python3 tools/package.py OUT_DIR [FOUNDATION_DIR]
       (foundation/ in the repository is included, and FOUNDATION_DIR adds more)

  opencif_evaluation.zip     the self-contained evaluation archive, which holds the package, tools, tests, presets,
                             examples, docs (with RESULTS, CONTACT_SHEET, OPTIONS, ALIGNMENT, LIMITATIONS, and
                             EVALUATION), the README, the requirements, every run's record and summary, the CFD-kit
                             tables, the build and test logs, the frozen provenance scripts, the foundation sources
                             (foundation/, and the files of FOUNDATION_DIR if it is given), and a MANIFEST of sha256
                             sums
  opencif_figures_N.zip      every run's proof figures (layout, zooms, physics) and CFD-kit windows, split by size
  opencif_masks_N.zip        every run's DXF photomask (and the array masks), split by size
  opencif_traces_N.zip       the micrograph-trace and whole-device figures, as PNG and PDF, split by size. They are
                             rc6's provenance/rc6_traces/figures/ and full_device/ and rc7's
                             provenance/rc7_traces/results/ with every folder below it, and their records, code, and
                             README stay in the evaluation archive. These figures show the papers' micrographs, so
                             they belong with the dissertation's evidence rather than in a public repository.

The Fig. 4 variants' builds in provenance/rc7_traces/variants/runs/ (in RC8, the rejected readings B, pathorder, and C,
start51, and the folded-head figure, fold20) are split three ways. Their DXFs go with the masks
(tests/test_trace_anchors.py layer 6 reads them), their figures with the figures, and their records in the evaluation
archive. Build pickles (*.pkl) are never packaged.
"""
import glob
import os
import sys
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIMIT_MB = 27.0                 # each zip under the 30 MiB upload limit


def add_tree(z, base, rel, skip=('__pycache__', '.pyc', '/runs/', '/logs/')):
    for dp, dn, fn in os.walk(os.path.join(base, rel)):
        for f in fn:
            p = os.path.join(dp, f)
            if any(s in p for s in skip):
                continue
            z.write(p, os.path.relpath(p, base))


def compressed_size(f):
    """Return an upper estimate of a file's size in the zip. PNG, PDF, and archive files count at their own size,
    because they are already compressed, and any other file is deflated here at level 6 (the zip uses 9, which is never
    larger), with the entry's headers added."""
    import zlib
    if f.lower().endswith(('.png', '.pdf', '.zip', '.jpg', '.jpeg', '.gz')):
        n = os.path.getsize(f)
    else:
        n = len(zlib.compress(open(f, 'rb').read(), 6))
    return n + 200 + 2 * len(os.path.relpath(f, ROOT))


def split_zip(out_dir, stem, files):
    """Write the files into stem_1.zip, stem_2.zip, and so on, each under LIMIT_MB. A zip is closed before a file that
    would take it over the limit, whereas rc5 checked after each file, so that one large figure could overshoot."""
    made, i, z, cur, size = [], 0, None, None, 0
    for f in sorted(files):
        est = compressed_size(f)
        assert est < LIMIT_MB * 1e6, (f, est)            # one file alone must fit
        if z is not None and size + est > LIMIT_MB * 1e6 * 0.98:
            z.close()
            z = None
        if z is None:
            i += 1
            cur = os.path.join(out_dir, '%s_%d.zip' % (stem, i))
            z = zipfile.ZipFile(cur, 'w', zipfile.ZIP_DEFLATED, compresslevel=9)
            made.append(cur)
            size = 22                                      # the end-of-archive record
        z.write(f, os.path.relpath(f, ROOT))
        z.fp.flush()
        # the archive so far, plus the directory
        size = os.path.getsize(cur) + 46 * len(z.infolist()) + sum(len(q.filename) for q in z.infolist())
    if z is not None:
        z.close()
    return made


def main(out_dir, foundation=None):
    import hashlib
    os.makedirs(out_dir, exist_ok=True)
    made = []
    p = os.path.join(out_dir, 'opencif_evaluation.zip')
    top = 'opencif_evaluation'
    entries = []                                          # (archive path, source path)
    trace_figs, var_masks, var_figs = [], [], []
    for rel in ('opencif', 'tools', 'tests', 'presets', 'examples', 'docs', 'reference_glossary', 'provenance', 'foundation'):
        for dp, dn, fn in os.walk(os.path.join(ROOT, rel)):
            for f in fn:
                q = os.path.join(dp, f)
                if '__pycache__' in q or q.endswith('.pyc'):
                    continue
                rq = os.path.relpath(q, ROOT)
                if f.endswith('.pkl'):
                    continue                      # the build's working arrays, which are not results
                if (rq.startswith(os.path.join('provenance', 'rc6_traces', 'figures') + os.sep) or
                        rq.startswith(os.path.join('provenance', 'rc6_traces', 'full_device') + os.sep) or
                        rq.startswith(os.path.join('provenance', 'rc7_traces', 'results') + os.sep)) and \
                        f.lower().endswith(('.png', '.pdf')):
                    trace_figs.append(q)          # the trace figures go to their own zips
                    continue
                if rq.startswith(os.path.join('provenance', 'rc7_traces', 'variants', 'runs') + os.sep) and \
                        f.lower().endswith(('.dxf', '.png')):
                    (var_masks if f.lower().endswith('.dxf') else var_figs).append(q)
                    continue
                entries.append((rq, q))
    for f in ('README.md', 'requirements.txt', 'requirements-tested.txt'):
        entries.append((f, os.path.join(ROOT, f)))
    entries.append(('EVALUATION.md', os.path.join(ROOT, 'docs', 'EVALUATION.md')))
    for f in sorted(glob.glob(os.path.join(ROOT, 'runs', '*', '*.json')) + glob.glob(os.path.join(ROOT, 'runs', '*', '*.md'))
                    + glob.glob(os.path.join(ROOT, 'runs', '*', '*.yaml'))
                    + glob.glob(os.path.join(ROOT, 'runs', '*', 'cfd_kit_*', '*.csv'))
                    + glob.glob(os.path.join(ROOT, 'runs', '*', 'cfd_kit_*', '*.md'))):
        if '/_t/' not in f:
            entries.append((os.path.relpath(f, ROOT), f))
    # The comparison logs that the documents cite and the rebuild log ship too (final document check, m6).
    for f in sorted(glob.glob(os.path.join(ROOT, 'logs', '*.log')) + glob.glob(os.path.join(ROOT, 'logs', 'batch_summary.txt'))
                    + glob.glob(os.path.join(ROOT, 'logs', 'compare_*.md')) + glob.glob(os.path.join(ROOT, 'logs', 'compare_*.json'))
                    + glob.glob(os.path.join(ROOT, 'logs', 'rebuild_rc*.txt'))):
        entries.append((os.path.relpath(f, ROOT), f))
    if foundation:
        for f in sorted(os.listdir(foundation)):
            entries.append(('foundation/' + f, os.path.join(foundation, f)))
    lines = []
    with zipfile.ZipFile(p, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for arc, src in entries:
            z.write(src, top + '/' + arc)
            lines.append('%s  %s' % (hashlib.sha256(open(src, 'rb').read()).hexdigest(), arc))
        z.writestr(top + '/MANIFEST.txt', 'sha256 of every file in this archive (paths relative to %s/)\n' % top
                   + '\n'.join(lines) + '\n')
    made.append(p)
    figs = [f for f in glob.glob(os.path.join(ROOT, 'runs', '*', '*.png')) + glob.glob(os.path.join(ROOT, 'runs', '*', 'cfd_kit_*', '*.png'))
            if '/_t/' not in f] + var_figs
    made += split_zip(out_dir, 'opencif_figures', figs)
    masks = [f for f in glob.glob(os.path.join(ROOT, 'runs', '*', '*.dxf')) + glob.glob(os.path.join(ROOT, 'runs', '*', 'cfd_kit_*', '*.dxf'))
             if '/_t/' not in f] + var_masks
    made += split_zip(out_dir, 'opencif_masks', masks)
    if trace_figs:
        made += split_zip(out_dir, 'opencif_traces', trace_figs)
    for m in made:
        print('%-60s %8.2f MB' % (os.path.basename(m), os.path.getsize(m) / 1e6))
        assert os.path.getsize(m) < LIMIT_MB * 1e6, m


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'dist'), sys.argv[2] if len(sys.argv) > 2 else None)
