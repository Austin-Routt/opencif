"""Write notebooks/opencif_traces.ipynb (unexecuted), the companion of the tutorial. It assembles the four device
packages from the repository and the publishers' figures, and builds, traces, measures, and renders each of them.
Usage: python make_traces_v1.py OUT.ipynb"""
import sys
import nbformat as nbf

cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s.strip('\n')))
code = lambda s: cells.append(nbf.v4.new_code_cell(s.strip('\n')))

md(r"""
# opencif 1.0.0 on the papers' micrographs
### The four device packages: photomasks, micrograph traces, and key figures

This notebook reads the publishers' figures, and its outputs embed their pixels. Gifford et al. (2014, Figs. 4 and 6)
and Dinh et al. (2024, Figs. 2 and 3) are © the Royal Society of Chemistry, and the articles carry no Creative Commons
license; their figures are reproduced with the Royal Society of Chemistry's permission. Xia et al. (2016) is published
under CC BY 4.0, and its attribution is in the Xia kit's `inputs/SOURCE.md`. The figures are in `notebooks/figures/`,
outside opencif's MIT License (`NOTICE.md`), and the notebook's outputs carry the same terms.

This notebook is the companion of `opencif_tutorial.ipynb`. For each published device, it assembles the device package
from this repository (the package, the presets, the trace code, and the kit, each checked by its sha256) and from the
figures in `notebooks/figures/`, and it runs the package's script (`devices/README.md`). The script runs five steps.
First, it checks the code stamp against the release's and every file of the package against its `MANIFEST.json`.
Second, it builds every photomask of the package and checks that each is byte-identical to the release's photomask,
which is RC8's. Third, it registers each photomask on the paper's micrograph, with the scale taken from the printed
bars or widths, the rotation measured on the image, and two independent registrations compared, and it checks that
every value of each trace record equals RC8's record. Fourth, it measures the alignment (widths, offsets, and null
ratios for Fig. 6 and for Xia, Fig. 3a's boxes for Xia, the taper-direction test V13 followed by the measurement for
Gifford's Fig. 4, and the null ratios of the trace step for Dinh). Finally, it renders the figures, which are each
micrograph with the photomask's edges over it (with no text added) and the whole device with the micrographs placed by
their registrations. Any departure from RC8's records stops a package with exit code 2.

The figures were extracted from the articles with poppler's `pdfimages`, which writes the publishers' embedded bitmaps
without resampling, and `RASTERS` (below) reads them from `notebooks/figures/`. To use figures extracted from your own
copies of the articles, set the environment variable `OPENCIF_PAPER_RASTERS` to their folder.

```
pdfimages -png -f 5 -l 5 "<Dinh 2024 PDF>"    dinh5    # dinh5-000.png  Fig. 2
pdfimages -png -f 6 -l 6 "<Dinh 2024 PDF>"    dinh6    # dinh6-000.png  Fig. 3
pdfimages -png -f 5 -l 5 "<Gifford 2014 PDF>" gif5     # gif5-001.png   Fig. 4 (gif5-000.png is Fig. 3, not used)
pdfimages -png -f 7 -l 7 "<Gifford 2014 PDF>" gif7     # gif7-000.png   Fig. 6
pdfimages -png -f 5 -l 5 "<Xia 2016 PDF>"     xia5     # xia5-000.png   Fig. 3
```

A device takes about 45 to 75 minutes, most of it in the builds, and the notebook runs two devices at a time, so the
whole run takes about 2.5 to 3 hours on two cores. A device whose package finished in an earlier run of this notebook
(its run log ends with a complete run) is not run again. Start Jupyter in the repository's `notebooks/` folder.
""")
code(r"""
import os, sys, json, time, hashlib, subprocess, glob
from pathlib import Path
from IPython.display import display, Markdown, Image
HERE = Path.cwd()
TREE = next(p for p in [HERE, *HERE.parents] if (p / 'opencif' / '__init__.py').exists())   # the repository
sys.path.insert(0, str(TREE))
from opencif import pipeline as PL
RELEASE = json.load(open(TREE / 'validation' / 'RELEASE.json'))
RASTERS = Path(os.environ.get('OPENCIF_PAPER_RASTERS', HERE / 'figures'))
OUT = HERE / 'traces_out'
OUT.mkdir(exist_ok=True)
DEVICES = ['dinh2024_v1', 'gifford2014_fig4', 'gifford2014_fig6', 'xia2016_rs7']


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


print('opencif %s, code stamp %s; the release\'s code stamp %s (%s)' % (
    PL.CODE['version'], PL.CODE['sha256'], RELEASE['code_stamp'],
    'equal' if PL.CODE['sha256'] == RELEASE['code_stamp'] else 'DIFFERENT, so the packages will stop at their check'))
need = {}
for dev in DEVICES:
    man = json.load(open(TREE / 'devices' / dev / 'MANIFEST.json'))
    for rel, src in man['files'].items():
        if 'raster' in src:
            need[Path(rel).name] = (src['raster'], dev)
print('figures folder:', RASTERS)
ok = True
for name, (h, dev) in sorted(need.items()):
    p = RASTERS / name
    got = sha256(p) if p.exists() else None
    ok &= got == h
    print('  %-14s for %-17s %s' % (name, dev, 'present, sha256 as expected' if got == h else
                                     ('MISSING' if got is None else 'present but different (%s...)' % got[:16])))
print('all figures present and verified:', ok)
""")
md(r"""
## Running the four packages

Each package is assembled and run in its own process by `python tools/device_package.py <name> --rasters ...
--out traces_out`, two at a time, with its output in `traces_out/<name>.log`. The package's own record of every step is
`traces_out/opencif_<name>/out/run_log.txt`.
""")
code(r"""
def finished(dev):
    '''Return True when the package's last run (from its last "== check") ended with every step done.'''
    log = OUT / ('opencif_' + dev) / 'out' / 'run_log.txt'
    if not log.exists():
        return False
    lines = log.read_text().splitlines()
    starts = [i for i, ln in enumerate(lines) if '== check' in ln]
    tail = lines[starts[-1]:] if starts else lines
    done = [ln for ln in tail if 'done:' in ln]
    return bool(done) and all(k in done[-1] for k in ('check', 'build', 'trace', 'render'))


EXPECTED_MIN = {'dinh2024_v1': 50, 'gifford2014_fig4': 70, 'gifford2014_fig6': 60, 'xia2016_rs7': 75}   # on two busy cores
BUDGET_MIN = float(os.environ.get('OPENCIF_TRACES_BUDGET_MIN', 'inf'))   # optional; no device starts that would end after it


def run_devices(devs, workers=2):
    todo = [d for d in devs if not finished(d)]
    for d in devs:
        if d not in todo:
            print('%-17s already complete in an earlier run (its run log ends with a complete run)' % d)
    running, t0 = {}, time.time()
    while todo or running:
        while todo and len(running) < workers:
            if (time.time() - t0) / 60 + EXPECTED_MIN[todo[0]] > BUDGET_MIN:
                print('%-17s not started, because the time budget (%g min) would be exceeded; run this cell again'
                      % (todo[0], BUDGET_MIN))
                todo.pop(0)
                continue
            d = todo.pop(0)
            log = open(OUT / (d + '.log'), 'a')
            p = subprocess.Popen([sys.executable, '-I', str(TREE / 'tools' / 'device_package.py'), d, '--rasters',
                                  str(RASTERS), '--out', str(OUT)], cwd=str(TREE), stdout=log, stderr=subprocess.STDOUT)
            running[d] = (p, log, time.time())
            print('%-17s started %s' % (d, time.strftime('%H:%M:%S')))
        for d, (p, log, ts) in list(running.items()):
            if p.poll() is not None:
                log.close()
                print('%-17s exit %d after %.0f min' % (d, p.returncode, (time.time() - ts) / 60))
                del running[d]
        time.sleep(15)
    print('finished after %.0f min' % ((time.time() - t0) / 60))


if ok:
    run_devices(DEVICES)
else:
    print("the figures in RASTERS are missing or differ from the kits' sha256")
""")
md(r"""
## What each package checked

The lines below come from the packages' own run logs. They show every input's sha256, each photomask against the
release's (which is RC8's), every trace record against RC8's record, and, for Fig. 4, the taper-direction test against
its recorded verdicts.
""")
code(r"""
KEYS = ('byte-identical', 'every value of', 'taper (', 'exit', 'done:', 'FAIL', 'Error', 'differs', 'package:')
for dev in DEVICES:
    log = OUT / ('opencif_' + dev) / 'out' / 'run_log.txt'
    display(Markdown('### %s' % dev))
    if not log.exists():
        print('not run')
        continue
    lines = log.read_text().splitlines()
    runs = [i for i, ln in enumerate(lines) if '== check' in ln]
    sel = [ln for ln in lines[runs[-1] if runs else 0:] if any(k in ln for k in KEYS)]
    print('\n'.join(ln[:230] for ln in sel))
""")
md(r"""
## The figures

The cell below previews each package's main figures, whose files at full resolution are in
`traces_out/opencif_<name>/out/`. They are the zoomed micrographs with the photomask's edges over them (cyan where the
registration makes a claim and amber where it does not) and the whole device with the micrographs placed by their
registrations, together with the taper-direction figure for Gifford's Fig. 4 (widths only).
""")
code(r"""
def preview(path, width=1100):
    '''Show a PNG scaled to `width` pixels; the file on disk keeps its full resolution.'''
    from PIL import Image as PILImage
    import io
    PILImage.MAX_IMAGE_PIXELS = None
    im = PILImage.open(path)
    if im.width > width:
        im = im.resize((width, max(1, round(im.height * width / im.width))), PILImage.LANCZOS)
    buf = io.BytesIO()
    im.convert('RGB').save(buf, format='PNG', optimize=True)
    display(Image(data=buf.getvalue()))
    print('(preview; full resolution: %s)' % os.path.relpath(path, HERE))


SHOW = {'dinh2024_v1': ['out/figure_aligned/fig_zoom_overlay.png', 'out/figure_aligned/fig_full_device_stacked.png'],
        'gifford2014_fig4': ['out/fig_zoom_overlay.png', 'out/fig_taper_direction.png', 'out/fig_full_device_stacked.png'],
        'gifford2014_fig6': ['out/fig_zoom_overlay.png', 'out/fig_full_device_stacked.png'],
        'xia2016_rs7': ['out/fig_zoom_overlay.png', 'out/fig_drawing_overlay.png']}
for dev, files in SHOW.items():
    for f in files:
        p = OUT / ('opencif_' + dev) / f
        if p.exists():
            display(Markdown('%s, `%s`' % (dev, f)))
            preview(p)
""")
md(r"""
## The registrations in numbers

Each package's `out/registration.json` gives, for each panel, the scale, the rotation, the agreement of the two
registrations (V1), the median and 90th-percentile distance from the photomask's edges to the image's outlines, the
coverage (the share of eligible samples on which an outline was found), and the null ratios, which measure how much
worse the fit becomes when the mask is shifted by half a row pitch or by 5 µm (a ratio below 1.3 is flagged). RC8's
anchors decide which panels carry a claim and which are only recorded (`provenance/rc7_traces/README.md`).
""")
code(r"""
anchors = json.load(open(TREE / 'provenance' / 'rc7_traces' / 'trace_anchors.json'))
for dev in DEVICES:
    reg = OUT / ('opencif_' + dev) / 'out' / 'registration.json'
    if not reg.exists():
        continue
    R = json.load(open(reg))
    display(Markdown('### %s' % dev))

    def rows(panels, label):
        for k, p in panels.items():
            m, nr = p.get('metric') or {}, p.get('null_ratio') or {}
            print('  %-24s %-8s median %.2f um, p90 %.2f um, coverage %s, null ratios metric %s / censored %s / chamfer %s'
                  % (label, k, m.get('median_abs_um', float('nan')), m.get('p90_abs_um', float('nan')),
                     ('%.3f' % m['coverage']) if m.get('coverage') is not None else '-',
                     *[('%.2f' % nr[q]) if nr.get(q) is not None else '-' for q in ('metric', 'metric_censored', 'chamfer')]))
    for key, val in R.items():
        if isinstance(val, dict) and isinstance(val.get('panels'), dict):
            rows(val['panels'], key)
    if isinstance(R.get('panels'), dict):
        rows(R['panels'], 'device')
print('RC8 trace anchors: %d panels, %d null-ratio flags, code stamp %s' % (len(anchors['panels']), len(anchors['null_flags']),
                                                                         anchors['code_stamp'].get('sha256')))
""")
md(r"""
## The record

The last cell lists every file that the packages wrote, with its sha256, in `traces_out/MANIFEST.txt`, and checks each
package's own `out/SHA256SUMS`.
""")
code(r"""
lines = ['opencif %s, code stamp %s (the release: %s)' % (PL.CODE['version'], PL.CODE['sha256'], RELEASE['code_stamp'])]
for dev in DEVICES:
    sums = OUT / ('opencif_' + dev) / 'out' / 'SHA256SUMS'
    if sums.exists():
        ok_ = subprocess.run(['sha256sum', '-c', '--quiet', 'SHA256SUMS'], cwd=sums.parent, capture_output=True, text=True)
        lines.append('%s: out/SHA256SUMS %s (%d files)' % (dev, 'verifies' if ok_.returncode == 0 else 'DOES NOT VERIFY',
                                                          len(sums.read_text().splitlines())))
files = sorted(p for p in OUT.rglob('*') if p.is_file() and p.name != 'MANIFEST.txt' and '.mplconfig' not in p.parts)
(OUT / 'MANIFEST.txt').write_text('\n'.join(lines + ['%s  %s' % (sha256(p), p.relative_to(OUT)) for p in files]) + '\n')
print('\n'.join(lines))
print('%d files listed in traces_out/MANIFEST.txt' % len(files))
""")

nb = nbf.v4.new_notebook()
nb['cells'] = cells
nb['metadata'] = {'kernelspec': {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'},
                  'language_info': {'name': 'python'}}
nbf.write(nb, sys.argv[1])
print('wrote %s: %d cells' % (sys.argv[1], len(cells)))
