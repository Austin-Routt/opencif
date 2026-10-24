"""devicekit.py: the steps that every opencif 1.0.0 device package runs, in the order in which it runs them.

  check   The package's copy of opencif is the release (its code stamp is the one in tree/validation/RELEASE.json), and
          every file of the package is the one that MANIFEST.json names (sha256), including the trace code and the
          publishers' figures.
  build   The photomask is built by opencif's pipeline.run on the device's preset (the generator's input), and the
          DXF's sha256 is compared with the release's photomask, which is RC8's byte for byte (reference/).
  trace   The micrograph registrations are made by the trace code (tree/provenance/rc7_traces/code) on the photomask
          just built, and every value of the new registration record is compared with RC8's trace record
          (reference/).
  render  Three figures are drawn with no text added, so that the micrographs keep the publishers' own labels. They
          are the zoomed micrograph panels with the mask's edges over them, the whole device with the micrographs
          placed by their registrations, and that drawing cut into consecutive rows stacked from top to bottom, each
          as a vector PDF (micrographs embedded at their native pixels and edges as vectors) and as a PNG.

Nothing here measures anything new. The build is the generator's, the registrations are the trace code's, and the
figures draw what those two produced. The references are RC8's own records (its rebuild and its traces), which were
made independently of this package and which opencif 1.0.0 reproduces, because its photomasks are RC8's byte for byte.
The device script (make_<device>.py) holds everything that is specific to a device.
"""
import hashlib
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
# matplotlib keeps its cache and configuration inside the package and outside out/, because nothing is written outside
# this folder and the cache is not part of the results (rc6 package check, finding 29). The variable is assigned rather
# than defaulted, since an inherited one may lie outside the folder.
os.environ['MPLCONFIGDIR'] = os.path.join(HERE, '.mplconfig')
os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)

import numpy as np   # noqa: E402

TREE = os.path.join(HERE, 'tree')
# The trace code stays at its place in the tree, because some of its scripts find the tree's data files relative to
# their own path.
TRACE = os.path.join(TREE, 'provenance', 'rc7_traces', 'code')
RELEASE_JSON = os.path.join(TREE, 'validation', 'RELEASE.json')
STAMP = json.load(open(RELEASE_JSON))['code_stamp'] if os.path.exists(RELEASE_JSON) else None   # the release's stamp
CYAN = '#00d5ff'            # the mask's edges over a micrograph (the trace color)
ORANGE = '#ff6a00'          # a placed micrograph's frame on the whole device
AMBER = '#ffb000'           # the mask's edges over a micrograph where the panel makes no claim (the no-claim color)
FLUID, ISLAND, EDGE = '#cfe3f7', '#ffffff', '#1b4f9c'
PDF_METADATA = {'CreationDate': None}   # no creation date in the PDFs, so a rerun gives byte-identical files

LOG = []


def say(msg):
    """Print one time-stamped line of the run's log and keep it in LOG, which the device script appends to
    out/run_log.txt."""
    line = '%s  %s' % (time.strftime('%Y-%m-%d %H:%M:%S'), msg)
    print(line, flush=True)
    LOG.append(line)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def yaml_load(fh):
    import yaml
    return yaml.safe_load(fh)


def fail(msg):
    """Log the failure and stop the run with exit code 2."""
    say('FAILED: ' + msg)
    raise SystemExit(2)


# ------------------------------------------------------------------------------------------------------------ check
def use_tree():
    """Put the package's copy of the release (tree/) first on the import path, together with the trace code and the
    measurement scripts of provenance/rc7_scripts that the device scripts import."""
    for p in (os.path.join(TREE, 'provenance', 'rc7_scripts'), TRACE, TREE):
        if p not in sys.path:
            sys.path.insert(0, p)


def check_packages():
    """Check that the Python packages that opencif needs are present, and stop if one is missing. A package at a
    version other than the tested one is noted in the log rather than refused."""
    import importlib
    have, missing = {}, []
    for mod, name in (('numpy', 'numpy'), ('scipy', 'scipy'), ('shapely', 'shapely'), ('matplotlib', 'matplotlib'),
                      ('yaml', 'pyyaml'), ('ezdxf', 'ezdxf'), ('PIL', 'pillow')):
        try:
            have[name] = importlib.import_module(mod).__version__
        except ImportError:
            missing.append(name)
    if missing:
        fail('python %s lacks %s: install tree/requirements-tested.txt (pip install -r tree/requirements-tested.txt)'
             % (sys.version.split()[0], ', '.join(missing)))
    say('python %s; %s' % (sys.version.split()[0], ', '.join('%s %s' % kv for kv in have.items())))
    tested = {}
    with open(os.path.join(TREE, 'requirements-tested.txt')) as fh:
        req = fh.read().splitlines()
    for line in req:
        line = line.split('#')[0].strip()
        if '==' in line:
            k, v = line.split('==')
            tested[k.strip().lower()] = v.strip()
    other = ['%s %s (tested %s)' % (k, have.get(k), v) for k, v in tested.items() if have.get(k) != v]
    if other:
        say('NOTE: not the tested versions: %s.  The photomask and the registrations are still checked against the '
            'release\'s records; the figures may differ in their bytes' % '; '.join(other))
    else:
        say('every package at the version in tree/requirements-tested.txt')


def check_code():
    """Check that the package's opencif is the release, because its code stamp equals the one in
    tree/validation/RELEASE.json, and that every file of the package has the sha256 that MANIFEST.json gives."""
    check_packages()
    use_tree()
    from opencif import pipeline
    if STAMP is None or pipeline.CODE['sha256'] != STAMP:
        fail('opencif code stamp %s, expected the release\'s %s (tree/validation/RELEASE.json)'
             % (pipeline.CODE['sha256'], STAMP))
    say('opencif %s, code stamp %s (%d hashed files): the release' % (pipeline.CODE['version'], pipeline.CODE['sha256'],
                                                                     pipeline.CODE['files']))
    # Every file that the package was assembled from is checked again here, so that a package edited after its
    # assembly is found. The trace code is among these files.
    with open(os.path.join(HERE, 'MANIFEST.json')) as fh:
        man = json.load(fh)
    bad, n, n_trace = [], 0, 0
    for rel, src in man['files'].items():
        want = src.get('sha256') or src.get('raster')
        n += 1
        n_trace += rel.startswith('tree/provenance/rc7_traces/code/')
        if sha256(os.path.join(HERE, *rel.split('/'))) != want:
            bad.append(rel)
    if bad:
        fail('files differ from MANIFEST.json: %s' % ', '.join(bad[:20]))
    say('package: %d files (%d of the trace code), every sha256 as MANIFEST.json lists it' % (n, n_trace))


def check_inputs(files):
    """Check the sha256 of each input file against its expected value, and stop at the first that differs."""
    for path, want in files.items():
        got = sha256(path)
        if got != want:
            fail('%s: sha256 %s, expected %s' % (os.path.relpath(path, HERE), got, want))
        say('input %s: sha256 %s (as expected)' % (os.path.relpath(path, HERE), got[:16]))


def check_preset(preset, inputs, skip=('name', 'description', 'paper')):
    """Check the generator input (the preset) against the device script's INPUTS table in both directions. Every value
    in the table must be the YAML's, with the same type, and every value of the YAML outside the free-text sections in
    `skip` must be in the table, so that the table is the whole input."""
    import yaml
    spec = yaml.safe_load(open(preset))

    def leaves(d, p=''):
        for k, v in d.items():
            q = '%s.%s' % (p, k) if p else str(k)
            if isinstance(v, dict):
                yield from leaves(v, q)
            else:
                yield q, v
    have = dict(leaves(spec))
    table = {}
    for path, value, code, src in inputs:
        if code not in 'PMDAG' or len(code) != 1:
            fail('INPUTS %s: code %r is not one of P M D A G' % (path, code))
        if path in table:
            fail('INPUTS lists %s twice' % path)
        table[path] = value
        if path not in have:
            fail('INPUTS %s is not in the preset' % path)
        if have[path] != value or type(have[path]) is not type(value):
            fail('input %s is %r in the preset, %r in INPUTS' % (path, have[path], value))
    missing = [q for q in have if q.split('.')[0] not in skip and q not in table]
    if missing:
        fail('the preset sets %s, which INPUTS does not list' % ', '.join(missing))
    say('generator input: the preset holds exactly the %d values of INPUTS (P %d, M %d, D %d, A %d, G %d; free text '
        'under %s aside)' % (len(inputs), *[sum(1 for i in inputs if i[2] == c) for c in 'PMDAG'], '/'.join(skip)))


# ------------------------------------------------------------------------------------------------------------ build
def build(preset, outdir, expect_dxf_sha256, expect_exit):
    """Build the preset with the generator and return (dxf, record path, record). The DXF must be byte-identical to
    the release's photomask of the same preset, which is RC8's (reference/<name>_dxf_sha256.txt), and the exit code
    must be the release's."""
    use_tree()
    from opencif import pipeline as PL
    import yaml
    name = yaml.safe_load(open(preset))['name']
    say('build: opencif pipeline.run(%s) -> %s' % (os.path.relpath(preset, HERE), os.path.relpath(outdir, HERE)))
    t0 = time.time()
    rec = PL.run(preset, outdir=outdir, verbose=True, render=False)
    a = rec['assessment']
    say('build: %s in %.0f s; exit code %d; N %d gap rows (%d in the meander head), F %.5f' % (
        name, time.time() - t0, a['exit_code'], rec['design']['N'], rec['design']['NM'], rec['design']['F']))
    say('build: verdict: %s' % a['verdict'])
    if a['exit_code'] != expect_exit:
        fail('exit code %d, expected %d (the release\'s record\'s)' % (a['exit_code'], expect_exit))
    dxf, record = os.path.join(outdir, name + '.dxf'), os.path.join(outdir, name + '.json')
    pkl = os.path.join(outdir, name + '.pkl')
    if os.path.exists(pkl):
        os.remove(pkl)                         # the build's working arrays (14 MB), not part of the photomask
    got = sha256(dxf)
    if got != expect_dxf_sha256:
        fail('photomask %s: sha256 %s, expected the release\'s (and RC8\'s) %s' % (os.path.basename(dxf), got,
                                                                                  expect_dxf_sha256))
    say('build: photomask %s, sha256 %s: byte-identical to the release\'s photomask, which is RC8\'s'
        % (os.path.basename(dxf), got))
    with open(record) as fh:
        return dxf, record, json.load(fh)


def check_build_physics(rec, expect_declared=(), f_tol=1e-4):
    """Check what the frozen physics says about the drawn device, from the build record's network solve of the drawn
    photomask, and return physics_block(rec). The delivered filtrate fraction must be the target to within f_tol
    (1e-4, or the record's own Q2 tolerance of 1e-3 where the device script says why), no row may leave the cutoff band
    without a declared reason, the declared departures must be exactly the expected ones, and the drawing must be one
    fluid component with nothing lost and no invalid ring."""
    v = rec['validation']
    rc = v['realised_cutoff']
    F_t, F_n = rec['input']['targets']['filtrate_fraction'], rc['F_net']
    if abs(F_n - F_t) > f_tol:
        fail('network F %.5f against the target %.5f (tolerance %g)' % (F_n, F_t, f_tol))
    if rc.get('undeclared_rows_outside', 0):
        fail('%d rows leave the cutoff band without a declared reason' % rc['undeclared_rows_outside'])
    dec = sorted(rc.get('declared_rows', {}))
    if dec != sorted(expect_declared):
        fail('declared departures %s, expected %s' % (dec, sorted(expect_declared)))
    d = v['drc']
    if d['n_components'] != 1 or d['n_lost_open'] or d['n_lost_solid'] or d['invalid_rings']:
        fail('drawing check: %d fluid components, %d / %d features lost, %d invalid rings' % (
            d['n_components'], d['n_lost_open'], d['n_lost_solid'], d['invalid_rings']))
    def _one(k):
        """Describe one declared departure in the log by its rows, its worst row, and its share of the filtrate."""
        o = rc['declared_by_option'].get(k) or {}
        fl = o.get('first_last') or []
        rows = ('rows %d-%d' % (fl[0], fl[1])) if len(fl) == 2 else 'no rows outside'
        wo = o.get('worst_over')
        return '%s (%s, %s, %.2f %% of the filtrate oversize)' % (
            k, rows, ('up to %.2f x D' % wo) if wo is not None else 'no row above D', 100 * (o.get('share_of_F_oversize') or 0.0))
    dtxt = ', '.join(_one(k) for k in dec)
    say('physics: network F %.5f (target %.5f); every row within 1 %% of D%s; one fluid component, nothing lost at '
        'the %.1f um drawing check' % (F_n, F_t, (' except the declared ' + dtxt) if dec else ' (no declared departure)',
                                       d['d_um']))
    ph = physics_block(rec)
    say('physics: %.2f %% of the filtrate through rows within 1 %% of D (%d rows outside); reverse flow %s; plateau '
        'u_min %.4f (%s); junction band: %s' % (
            100 * ph['share_of_F_within_1pc'], ph['rows_outside_1pc'],
            ('%.2f %% of the inlet flow over %d rows' % (100 * ph['reverse_flow_share_of_inlet'], ph['reverse_flow_rows']))
            if ph['reverse_flow_rows'] else 'none', ph['plateau_u_min'], 'inside' if ph['plateau_ok'] else 'below 0.035',
            ph['junction_band_text']))
    return ph


def physics_block(rec):
    """Return the build record's physics of the drawn device for alignment.json and the README. The block holds the
    share of F within 1% of D, the rows outside, the reverse flow, the plateau, the cutoff band below it, the declared
    departures, and the model band's junction terms (rc6 package check, findings 11, 12, and 16)."""
    v = rec['validation']
    rc, g, ms = v['realised_cutoff'], v['gates'], v.get('model_scenarios') or {}
    jb = [dict(ms[k], end=k[len('junction_band_'):]) for k in ('junction_band_lo', 'junction_band_hi') if isinstance(ms.get(k), dict)]
    rf = rc.get('reverse_flow_rows') or []
    out = dict(
        F_net=rc['F_net'], F_target=rec['input']['targets']['filtrate_fraction'], D_target_um=rc['D_target_um'],
        share_of_F_within_1pc=rc['share_of_F_within_1%'], rows_outside_1pc=rc['rows_outside_1%'],
        share_of_F_oversize_1pc=rc['share_of_F_oversize_1%'], D_realised_max_um=rc.get('D_realised_max_um'),
        row1_cutoff_over_D=rc.get('row1_cutoff_over_D'),
        reverse_flow_rows=len(rf), reverse_flow_share_of_inlet=rc.get('reverse_flow_share_of_inlet') or 0.0,
        declared_by_option=rc.get('declared_by_option') or {}, undeclared_rows_outside=rc.get('undeclared_rows_outside', 0),
        plateau_u_min=g.get('plateau_u_min'), plateau_ok=g.get('plateau_ok'),
        cutoff_band_um=(g.get('cutoff_band_um') or {}).get('band_um') if not g.get('plateau_ok') else None,
        junction_band=[dict((k, b.get(k)) for k in ('end', 'F', 'rows_newly_outside_tol', 'rows_losing_their_cutoff',
                                                     'share_of_F_in_those_rows', 'max_cutoff_change', 'max_change_row'))
                       for b in jb],
        verdict=rec['assessment']['verdict'], exit_code=rec['assessment']['exit_code'],
        N=rec['design']['N'], NM=rec['design']['NM'],
        footprint_mm=(v.get('drawn_parameters') or {}).get('footprint_mm'))
    out['junction_band_text'] = '; '.join(
        'at the terms\' %s end %d rows newly outside the 1 %% band, %d losing their cutoff, %.2f %% of F, max change '
        '%.3f of D (row %s)' % (
            {'lo': 'low', 'hi': 'high'}.get(b['end'], b['end']), b['rows_newly_outside_tol'] or 0, b['rows_losing_their_cutoff'] or 0, 100 * (b['share_of_F_in_those_rows'] or 0.0),
            b['max_cutoff_change'] or 0.0, b['max_change_row']) for b in out['junction_band']) or 'not computed'
    return out


# ------------------------------------------------------------------------------------------------------------ trace
SKIP_KEYS = {'raster', 'mask'}             # file paths, which differ between machines (skipped only as strings)


def compare(new, ref, path='', diffs=None, rel=1e-9, ab=1e-9, extra=None):
    """Compare a new registration record with the reference record and return the list of differences. Every leaf of
    the reference must be in the new record and equal to it (numbers to within 1e-9 absolute plus 1e-9 relative, and
    strings, booleans, and None exactly), and a field that only the new record has goes to `extra`. A key in SKIP_KEYS
    is skipped only when its value is a string (a file path), so a number under such a key is still compared (rc6
    package check, finding 9)."""
    diffs = [] if diffs is None else diffs
    extra = [] if extra is None else extra
    if isinstance(ref, dict):
        if not isinstance(new, dict):
            diffs.append(path + ': not a mapping')
            return diffs
        for k in sorted(set(ref) | set(new), key=str):
            if k in SKIP_KEYS and isinstance(ref.get(k, new.get(k)), str):
                continue
            if k not in ref:
                extra.append('%s/%s = %r' % (path, k, new[k]))
                continue
            if k not in new:
                diffs.append('%s/%s: missing from the new record' % (path, k))
                continue
            compare(new[k], ref[k], path + '/' + str(k), diffs, rel, ab, extra)
    elif isinstance(ref, list):
        if not isinstance(new, list) or len(new) != len(ref):
            diffs.append(path + ': list length')
            return diffs
        for i, (a, b) in enumerate(zip(new, ref)):
            compare(a, b, '%s[%d]' % (path, i), diffs, rel, ab, extra)
    elif isinstance(ref, bool) or ref is None or isinstance(ref, str):
        if new != ref:
            diffs.append('%s: %r against the reference %r' % (path, new, ref))
    elif isinstance(ref, (int, float)):
        if not isinstance(new, (int, float)) or (math.isnan(ref) != math.isnan(new)) or (
                not math.isnan(ref) and abs(new - ref) > ab + rel * abs(ref)):
            diffs.append('%s: %r against the reference %r' % (path, new, ref))
    return diffs


def jsonable(o):
    """Convert a numpy scalar or array to the Python value that json can write."""
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    raise TypeError(type(o))


def check_trace(res, ref_path, out_path):
    """Write the new registration record to out_path, compare every value of it with RC8's trace record (reference/),
    and stop the run if any value differs."""
    res = json.loads(json.dumps(res, default=jsonable))
    with open(out_path, 'w') as fh:
        json.dump(res, fh, indent=1)
    with open(ref_path) as fh:
        ref = json.load(fh)
    extra = []
    diffs = compare(res, ref, extra=extra)
    if diffs:
        for d in diffs[:30]:
            say('  differs: ' + d)
        fail('%d values of the registration record differ from RC8\'s %s' % (len(diffs), os.path.basename(ref_path)))
    n = len(json.dumps(ref))
    say('trace: every value of RC8\'s %s is in the new registration record and equal (%d characters of '
        'record compared; file paths aside)' % (os.path.basename(ref_path), n))
    for e in extra:
        say('trace: recorded now, not in the reference record: %s' % e)
    return res


# ------------------------------------------------------------------------------------------------------------ geometry
class Placement:
    """A panel's stored registration (map_chamfer), which maps a point (x, y) of the mask in µm to the image's pixels
    by image px = origin + x ex + y ey."""

    def __init__(self, rec):
        m = rec['map_chamfer']
        self.o = np.array(m['origin_px'], float)
        self.M = np.c_[np.array(m['ex_px_per_um'], float), np.array(m['ey_px_per_um'], float)]
        self.B = np.linalg.inv(self.M)

    def to_px(self, xy):
        return self.o + (self.M @ np.asarray(xy, float).T).T

    def to_mask(self, px):
        return (self.B @ (np.asarray(px, float) - self.o).T).T

    def affine_px_to_mask(self):
        from matplotlib.transforms import Affine2D
        t = -self.B @ self.o
        return Affine2D(np.array([[self.B[0, 0], self.B[0, 1], t[0]], [self.B[1, 0], self.B[1, 1], t[1]], [0, 0, 1]]))

    def um_per_px(self):
        return 1.0 / math.sqrt(abs(np.linalg.det(self.M)))

    def device_axis_in_image_deg(self):
        """Return the direction of the mask's x axis (along the device) in the image, in degrees from the image's
        +column axis toward its +row axis. The rows run down the image, so 90 is straight down."""
        return math.degrees(math.atan2(self.M[1, 0], self.M[0, 0]))


def frame_of(rec):
    return [int(v) for v in (rec.get('micrograph_box') or rec['crop'])]


def edges_in_frame_px(mask, pl, frame, pad_um=20.0, claim=None):
    """Return the mask's edges inside a micrograph frame (c0, r0, c1, r1), in image px and clipped to the frame, as two
    lists, the edges that the panel's claim covers and the others. The function claim(edge) -> bool takes an edge in
    mask um, and claim=None claims every edge."""
    from shapely.geometry import LineString, box
    c0, r0, c1, r1 = frame
    corners = pl.to_mask(np.array([[c0, r0], [c1, r0], [c1, r1], [c0, r1]], float) - 0.5)
    x0, y0 = corners.min(axis=0) - pad_um
    x1, y1 = corners.max(axis=0) + pad_um
    clip = box(c0 - 0.5, r0 - 0.5, c1 - 0.5, r1 - 0.5)
    out, other = [], []
    for e in mask.edges(x0, x1, y0, y1):
        dest = out if claim is None or claim(e) else other
        g = LineString(pl.to_px(np.asarray(e))).intersection(clip)
        for part in (g.geoms if hasattr(g, 'geoms') else [g]):
            if not part.is_empty and part.geom_type == 'LineString':
                dest.append(np.asarray(part.coords))
    return out, other


def load_rgb(path):
    from PIL import Image
    return np.asarray(Image.open(path).convert('RGB'))


# ------------------------------------------------------------------------------------------------------------ render
def render_zoom(panels, mask, out_stem, magnify=8, gap_px=4, line_px=0.25, overlay=True):
    """Draw the panels side by side, each micrograph at its own pixels (neither resampled nor rotated) and magnified
    `magnify` times, with the mask's edges over it in cyan (amber where the panel makes no claim), `line_px` image
    pixels wide, and with no text. The call writes <stem>.png (nearest-neighbor magnification) and <stem>.pdf (the
    micrographs embedded at their native pixels and the edges as vectors)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    info = []
    widths = [p['frame'][2] - p['frame'][0] for p in panels]
    heights = [p['frame'][3] - p['frame'][1] for p in panels]
    W = sum(widths) + gap_px * (len(panels) - 1)
    H = max(heights)
    dpi = 100.0
    for fmt in ('png', 'pdf'):
        fig = plt.figure(figsize=(W * magnify / dpi, H * magnify / dpi), dpi=dpi)
        fig.patch.set_facecolor('white')
        x = 0
        for p, w, h in zip(panels, widths, heights):
            c0, r0, c1, r1 = p['frame']
            ax = fig.add_axes([x / W, 1.0 - h / H, w / W, h / H])
            ax.set_axis_off()
            img = p['rgb'][r0:r1, c0:c1]
            ax.imshow(img, interpolation='nearest' if fmt == 'png' else 'none',
                      extent=(c0 - 0.5, c1 - 0.5, r1 - 0.5, r0 - 0.5), zorder=1)
            if overlay:
                ax.add_collection(LineCollection(p['edges_px'], colors=CYAN, linewidths=line_px * magnify * 72.0 / dpi,
                                                 zorder=2, capstyle='round', joinstyle='round'))
                if p['edges_px_other']:
                    ax.add_collection(LineCollection(p['edges_px_other'], colors=AMBER,
                                                     linewidths=line_px * magnify * 72.0 / dpi, zorder=2,
                                                     capstyle='round', joinstyle='round'))
            ax.set_xlim(c0 - 0.5, c1 - 0.5)
            ax.set_ylim(r1 - 0.5, r0 - 0.5)
            ax.set_aspect('auto')              # the axes box already has the panel's aspect (one image px = magnify px)
            x += w + gap_px
        path = '%s.%s' % (out_stem, fmt)
        fig.savefig(path, dpi=dpi, facecolor='white', metadata=PDF_METADATA if fmt == 'pdf' else None)
        plt.close(fig)
        info.append(path)
    return dict(files=[os.path.basename(f) for f in info], magnify=magnify, width_px=W * magnify, height_px=H * magnify,
                line_width_image_px=line_px if overlay else None, panel_order=[p['key'] for p in panels])


def device_bbox(mask):
    """Return the bounding box of every polygon of the photomask as (x0, x1, y0, y1) in mask um."""
    a = np.array([[p[2], p[3], p[4], p[5]] for p in mask.polys])
    return float(a[:, 0].min()), float(a[:, 1].max()), float(a[:, 2].min()), float(a[:, 3].max())


def _draw_full(ax, mask, panels, view, lw_edge, lw_trace, lw_frame, interp, clip_pad_um=None):
    """Draw the photomask's polygons inside view, the placed micrographs with their frames, and the mask's edges over
    them. With clip_pad_um, the fluid outline is cut to the view plus this pad (outside the axes' own clip) before it
    is drawn, so that a view of part of the device does not carry the whole outline. render_stack uses the pad,
    whereas render_full draws the outline whole."""
    from matplotlib.collections import LineCollection, PolyCollection
    from matplotlib.patches import Polygon as MplPolygon
    x0, x1, y0, y1 = view
    ch, isl = [], []
    for (L, p, bx0, bx1, by0, by1) in mask.polys:
        if bx1 < x0 or bx0 > x1 or by1 < y0 or by0 > y1:
            continue
        if L == 'CHANNELS' and clip_pad_um is not None:
            import shapely
            from shapely.geometry import Polygon
            g = shapely.clip_by_rect(Polygon(p), x0 - clip_pad_um, y0 - clip_pad_um, x1 + clip_pad_um, y1 + clip_pad_um)
            for part in (g.geoms if hasattr(g, 'geoms') else [g]):
                if part.geom_type != 'Polygon' or part.is_empty:
                    continue
                if len(part.interiors):
                    fail('a clipped fluid outline has holes; the stacked figure cannot fill it')
                ch.append(np.asarray(part.exterior.coords))
            continue
        (ch if L == 'CHANNELS' else isl).append(p)
    ax.add_collection(PolyCollection(ch, facecolors=FLUID, edgecolors='none', zorder=0))
    ax.add_collection(PolyCollection(isl, facecolors=ISLAND, edgecolors='none', zorder=1))
    ax.add_collection(LineCollection([np.asarray(e) for e in mask.edges(x0, x1, y0, y1)], colors=EDGE,
                                     linewidths=lw_edge, zorder=2))
    for p in panels:
        c0, r0, c1, r1 = p['frame']
        cm = p['corners_um']
        if cm[:, 0].max() < x0 or cm[:, 0].min() > x1 or cm[:, 1].max() < y0 or cm[:, 1].min() > y1:
            continue
        img = p['rgb'][r0:r1, c0:c1]
        if interp != 'none' and img.shape[2] == 3:
            # The raster path resamples the affine-placed image into the axes. The image is therefore given an opaque
            # alpha channel, so that the pixels outside a rotated micrograph stay transparent, because an RGB array is
            # filled black there.
            img = np.dstack([img, np.full(img.shape[:2], 255, img.dtype)])
        im = ax.imshow(img, interpolation=interp, extent=(c0 - 0.5, c1 - 0.5, r1 - 0.5, r0 - 0.5), zorder=4)
        im.set_transform(p['pl'].affine_px_to_mask() + ax.transData)
        frame = MplPolygon(cm, closed=True, fill=False, edgecolor=ORANGE, linewidth=lw_frame, zorder=6)
        ax.add_patch(frame)
        lc = LineCollection([p['pl'].to_mask(e) for e in p['edges_px']], colors=CYAN, linewidths=lw_trace, zorder=5)
        ax.add_collection(lc)
        if p['edges_px_other']:
            ax.add_collection(LineCollection([p['pl'].to_mask(e) for e in p['edges_px_other']], colors=AMBER,
                                             linewidths=lw_trace, zorder=5))


def render_full(mask, panels, out_stem, pdf_magnification=32.0, png_um_per_px=2.0, tile_mm=40.0, margin_um=150.0):
    """Draw the whole device, every polygon of the photomask, with each micrograph placed by the inverse of its stored
    registration (nothing is re-fitted) and the mask's edges over it, and with no text. The call writes two files:
      <stem>.pdf  a vector page that shows the device at `pdf_magnification` times its size, with the micrographs
                  embedded at their native pixels (matplotlib draws an affine-placed image unsampled in a PDF with
                  interpolation 'none')
      <stem>.png  one image at `png_um_per_px` µm per pixel, rendered in tiles of `tile_mm` and then joined"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from PIL import Image
    X0, X1, Y0, Y1 = device_bbox(mask)
    X0, X1, Y0, Y1 = X0 - margin_um, X1 + margin_um, Y0 - margin_um, Y1 + margin_um
    Lx, Ly = X1 - X0, Y1 - Y0
    out = {}
    # ---- PDF (vector)
    um_per_in = 25400.0 / pdf_magnification
    fig = plt.figure(figsize=(Lx / um_per_in, Ly / um_per_in))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    pt_per_um = 72.0 / um_per_in
    _draw_full(ax, mask, panels, (X0, X1, Y0, Y1), lw_edge=0.5 * pt_per_um, lw_trace=0.3 * pt_per_um,
               lw_frame=2.0 * pt_per_um, interp='none')
    ax.set_xlim(X0, X1)
    ax.set_ylim(Y0, Y1)
    ax.set_aspect('equal')
    fig.savefig(out_stem + '.pdf', facecolor='white', metadata=PDF_METADATA)
    plt.close(fig)
    out['pdf'] = dict(file=os.path.basename(out_stem + '.pdf'), page_in=[round(Lx / um_per_in, 3), round(Ly / um_per_in, 3)],
                      magnification=pdf_magnification, view_um=[X0, X1, Y0, Y1])
    # ---- PNG (tiles joined)
    s = png_um_per_px
    Wpx, Hpx = int(math.ceil(Lx / s)), int(math.ceil(Ly / s))
    tile_px = int(tile_mm * 1000.0 / s)
    canvas = Image.new('RGB', (Wpx, Hpx), 'white')
    dpi = 100.0
    for t0 in range(0, Wpx, tile_px):
        tw = min(tile_px, Wpx - t0)
        xa, xb = X0 + t0 * s, X0 + (t0 + tw) * s
        fig = plt.figure(figsize=(tw / dpi, Hpx / dpi), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_axis_off()
        _draw_full(ax, mask, panels, (xa, xb, Y0, Y1), lw_edge=0.8 * 72.0 / dpi, lw_trace=1.0 * 72.0 / dpi,
                   lw_frame=2.0 * 72.0 / dpi, interp='antialiased')          # line widths in output pixels
        ax.set_xlim(xa, xb)
        ax.set_ylim(Y0, Y0 + Hpx * s)
        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
        plt.close(fig)
        canvas.paste(Image.fromarray(buf[:Hpx, :tw].copy()), (t0, 0))
    canvas.save(out_stem + '.png', optimize=False, compress_level=6)
    out['png'] = dict(file=os.path.basename(out_stem + '.png'), size_px=[Wpx, Hpx], um_per_px=s, view_um=[X0, X0 + Wpx * s, Y0, Y0 + Hpx * s])
    return out


def stack_cuts(Wpx, s, X0, panels, rows, pad_um=100.0):
    """Return where the stacked figure cuts the device, in columns of the full-device PNG (the PDF cuts at the same
    x). The device is cut into `rows` pieces of nearly equal length along x, and a cut that falls on a micrograph's
    frame (widened by pad_um on each side) moves to the frame's nearer side, so that no micrograph is split."""
    spans = sorted([(p['corners_um'][:, 0].min() - pad_um - X0) / s, (p['corners_um'][:, 0].max() + pad_um - X0) / s]
                   for p in panels)
    merged = []
    for a, b in spans:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    cuts = [0]
    for k in range(1, rows):
        c = k * Wpx / rows
        for a, b in merged:
            if a < c < b:
                c = a if c - a <= b - c else b
        cuts.append(int(round(c)))
    cuts.append(Wpx)
    if any(b <= a for a, b in zip(cuts, cuts[1:])):
        fail('stacked figure: the cuts %s are not increasing (too many rows for the frames)' % cuts)
    return cuts


def stack_layout(mask, panels, full, rows, margin_um=150.0, gap_um=400.0):
    """Return the stacked figure's geometry in the full-device PNG's pixels (the PDF uses the same numbers in um). The
    device is cut into `rows` consecutive pieces along x (stack_cuts), each piece is trimmed across to its own content
    (the fluid outline and any micrograph frame in it, with margin_um on each side), and the pieces are stacked with
    gap_um between them."""
    from shapely import clip_by_rect
    from shapely.geometry import Polygon
    s = full['png']['um_per_px']
    Wpx, Hpx = full['png']['size_px']
    X0, _, Y0, _ = full['png']['view_um']
    Yt = Y0 + Hpx * s
    cuts = stack_cuts(Wpx, s, X0, panels, rows)
    fluid = [Polygon(p) for L, p, *_ in mask.polys if L == 'CHANNELS']
    out = []
    for a, b in zip(cuts, cuts[1:]):
        xa, xb = X0 + a * s, X0 + b * s
        lo, hi = [], []
        for g in fluid:
            c = clip_by_rect(g, xa, Y0, xb, Yt)
            if not c.is_empty:
                lo.append(c.bounds[1])
                hi.append(c.bounds[3])
        for p in panels:
            cx = p['corners_um'][:, 0]
            if cx.max() > xa and cx.min() < xb:
                lo.append(p['corners_um'][:, 1].min())
                hi.append(p['corners_um'][:, 1].max())
        if not lo:
            fail('stacked figure: nothing to draw between x %.0f and %.0f um' % (xa, xb))
        r0 = max(0, int(math.floor((Yt - (max(hi) + margin_um)) / s)))      # PNG rows, top first
        r1 = min(Hpx, int(math.ceil((Yt - (min(lo) - margin_um)) / s)))
        out.append((a, b, r0, r1))
    gap_px = int(round(gap_um / s))
    W = max(b - a for a, b, _, _ in out)
    H = sum(r1 - r0 for _, _, r0, r1 in out) + (rows - 1) * gap_px
    return dict(rows=out, gap_px=gap_px, size_px=[W, H], aspect=W / H)


def render_stack(mask, panels, out_stem, full, rows, pdf_magnification=32.0):
    """Draw the whole device cut into `rows` consecutive pieces along the chip's x axis and stacked from top to bottom,
    starting from the left end. For a linear device this is the flow order, whereas a folded device shows every leg in
    each row. The rows abut exactly, with no overlap and nothing left out, and every row is drawn at one scale,
    aligned on the left, and trimmed across to its own content (stack_layout). The drawing is render_full's (`full`
    is its return value), and the call writes two files:
      <stem>.png  the full-device PNG's own pixels, cut and stacked (2 um/px, with nothing redrawn)
      <stem>.pdf  a vector page at `pdf_magnification` times the device, with the micrographs embedded at their
                  native pixels"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from PIL import Image
    s = full['png']['um_per_px']
    Wpx, Hpx = full['png']['size_px']
    X0, _, Y0, _ = full['png']['view_um']
    Yt = Y0 + Hpx * s
    lay = stack_layout(mask, panels, full, rows)
    Wst, Hst = lay['size_px']
    gap_px = lay['gap_px']
    # ---- PNG: the full-device PNG cut and stacked
    limit, Image.MAX_IMAGE_PIXELS = Image.MAX_IMAGE_PIXELS, None       # our own file, larger than PIL's default
    try:
        src = Image.open(os.path.join(os.path.dirname(out_stem), full['png']['file']))
        src.load()
    finally:
        Image.MAX_IMAGE_PIXELS = limit
    if src.size != (Wpx, Hpx):
        fail('stacked figure: the full-device PNG is %s, expected %s' % (src.size, (Wpx, Hpx)))
    canvas = Image.new('RGB', (Wst, Hst), 'white')
    top = 0
    for a, b, r0, r1 in lay['rows']:
        canvas.paste(src.crop((a, r0, b, r1)), (0, top))
        top += r1 - r0 + gap_px
    canvas.save(out_stem + '.png', optimize=False, compress_level=6)
    # ---- PDF: each row in its own axes, with the same view as the PNG's row
    um_per_in = 25400.0 / pdf_magnification
    pt_per_um = 72.0 / um_per_in
    Wum, Hum = Wst * s, Hst * s
    fig = plt.figure(figsize=(Wum / um_per_in, Hum / um_per_in))
    top = 0
    rows_out = []
    for a, b, r0, r1 in lay['rows']:
        xa, xb, ya, yb = X0 + a * s, X0 + b * s, Yt - r1 * s, Yt - r0 * s
        ax = fig.add_axes([0.0, 1.0 - (top + r1 - r0) / Hst, (b - a) / Wst, (r1 - r0) / Hst])
        ax.set_axis_off()
        _draw_full(ax, mask, panels, (xa, xb, ya, yb), lw_edge=0.5 * pt_per_um, lw_trace=0.3 * pt_per_um,
                   lw_frame=2.0 * pt_per_um, interp='none', clip_pad_um=50.0)
        ax.set_xlim(xa, xb)
        ax.set_ylim(ya, yb)
        ax.set_aspect('auto')              # the axes box already has the row's aspect (one um the same on every row)
        rows_out.append(dict(x_um=[round(xa, 3), round(xb, 3)], y_um=[round(ya, 3), round(yb, 3)],
                             png_box_px=[0, top, b - a, top + r1 - r0], full_png_box_px=[a, r0, b, r1]))
        top += r1 - r0 + gap_px
    fig.savefig(out_stem + '.pdf', facecolor='white', metadata=PDF_METADATA)
    plt.close(fig)
    where = {p['key']: next(k + 1 for k, r in enumerate(rows_out)
                            if r['x_um'][0] <= p['corners_um'][:, 0].min() and p['corners_um'][:, 0].max() <= r['x_um'][1])
             for p in panels}
    return dict(files=[os.path.basename(out_stem) + e for e in ('.pdf', '.png')], rows=rows, row=rows_out,
                row_length_mm=[round((r['x_um'][1] - r['x_um'][0]) / 1000.0, 3) for r in rows_out],
                row_height_mm=[round((r['y_um'][1] - r['y_um'][0]) / 1000.0, 3) for r in rows_out],
                gap_um=gap_px * s, panel_row=where,
                pdf=dict(page_in=[round(Wum / um_per_in, 3), round(Hum / um_per_in, 3)], magnification=pdf_magnification),
                png=dict(size_px=[Wst, Hst], um_per_px=s))


def prepare_panels(order, results, rasters, mask, claim=None):
    """Return, for each panel in order, its registration, its frame, the micrograph's pixels, and the mask's edges
    over the frame. The argument claim is a dict that maps a panel's key to claim(edge in mask um) -> bool for the
    panels whose claim covers only some edges, and the edges outside a claim are drawn in amber."""
    panels = []
    for key in order:
        rec = results[key]
        pl = Placement(rec)
        frame = frame_of(rec)
        c0, r0, c1, r1 = frame
        corners_um = pl.to_mask(np.array([[c0, r0], [c1, r0], [c1, r1], [c0, r1]], float) - 0.5)
        ec, eo = edges_in_frame_px(mask, pl, frame, claim=(claim or {}).get(key))
        panels.append(dict(key=key, name=rec['name'], rec=rec, pl=pl, frame=frame, rgb=rasters[rec['raster_sha256']],
                           edges_px=ec, edges_px_other=eo, corners_um=corners_um))
    return panels


NULL_RATIO_CRITERION = 1.3     # every panel's null ratios below it are flagged (rc6 package check, finding 7)


def null_ratio(r):
    """Return the registration's discrimination, computed as provenance/rc7_traces/code/make_anchors.py computes it.
    Each ratio is the smaller of the metric's medians at the two shifted registrations (half a pitch along, 5 um
    across) over its median at the registration, and the same is done for the chamfer score."""
    c = r['V2_chamfer']
    out = dict(chamfer=min(min(c['half_pitch_along_px']), min(c['five_um_across_px'])) / c['best_px'])
    mn = r.get('metric_null') or {}
    if mn:
        out['metric'] = min(mn['half_pitch_along']['median_abs_um'], mn['5um_across']['median_abs_um']) / r['metric']['median_abs_um']
        cr = mn.get('censored_ratios') or {}
        # The censored ratio judges the null over every eligible sample, and a sample with no outline is counted at the
        # window (`trace_lib.summarise`; A7).
        if cr:
            out['metric_censored'] = min(cr['half_pitch_along'], cr['5um_across'])
    return dict((k, round(v, 4)) for k, v in out.items())     # the coverage is in registration_summary's metric


def registration_summary(panels):
    """Return each panel's registration as out/registration.json records it, with its scale, its rotation, its map,
    its metric and null ratios, its chamfer scores, and the share of edge samples on which an outline was found."""
    rows = {}
    for p in panels:
        r = p['rec']
        rows[p['key']] = dict(
            name=r['name'], raster_sha256=r['raster_sha256'], frame_px_c0_r0_c1_r1=p['frame'],
            um_per_px=round(p['pl'].um_per_px(), 6),
            device_axis_in_image_deg=round(p['pl'].device_axis_in_image_deg(), 4),
            rotation_measured_on_image_deg=next((r.get('measured', {})[k] for k in ('rotation_deg', 'leg_direction_deg', 'rows_direction_deg')
                                                 if k in r.get('measured', {})), None),
            rotation_measured_key=next((k for k in ('rotation_deg', 'leg_direction_deg', 'rows_direction_deg') if k in r.get('measured', {})), None),
            map_chamfer=r['map_chamfer'], frame_corners_um=np.round(p['corners_um'], 3).tolist(),
            V1_feature_vs_chamfer_px=r['V1_feature_vs_chamfer_px'],
            metric=dict((k, r['metric'][k]) for k in ('median_abs_um', 'p90_abs_um', 'mean_signed_um', 'n', 'coverage', 'n_eligible')
                        if k in r['metric']),
            metric_null=dict((k, r['metric_null'][k]['median_abs_um']) for k in ('half_pitch_along', '5um_across')),
            null_ratio=null_ratio(r), null_ratio_criterion=NULL_RATIO_CRITERION,
            below_criterion=sorted(k for k, v in null_ratio(r).items() if v < NULL_RATIO_CRITERION),
            chamfer=dict(best_px=r['V2_chamfer']['best_px'], half_pitch_px=r['V2_chamfer']['half_pitch_along_px'],
                         five_um_across_px=r['V2_chamfer']['five_um_across_px']),
            V5_scale_scan=dict((k, r['V5_scale_scan'][k]) for k in ('best_along', 'best_across')) if 'V5_scale_scan' in r else None,
            outline_found_share_V3=round(sum(c['n'] for c in r['V3_by_class'].values()) /
                                         max(1, sum(c['total'] for c in r['V3_by_class'].values())), 4),
            outline_found_share_note='V3: edge samples with an image outline within the window / all classified edge '
                                     'samples (a lower bound: the denominator includes samples outside the crop or '
                                     'inside the exclusions)',
            edges_drawn=len(p['edges_px']), edges_drawn_no_claim=len(p['edges_px_other']))
    return rows


def write_sums(outdir, skip_dirs=('trace_work', 'verification')):
    """Write SHA256SUMS over out/, leaving out the trace code's working folder (raster crops, removed after every
    trace) and out/verification/, which is written after the clean-room run and is about that run."""
    lines = []
    for root, dirs, files in os.walk(outdir):
        dirs[:] = [x for x in dirs if not (root == outdir and x in skip_dirs)]
        for f in sorted(files):
            if f == 'SHA256SUMS':
                continue
            p = os.path.join(root, f)
            lines.append('%s  %s' % (sha256(p), os.path.relpath(p, outdir)))
    open(os.path.join(outdir, 'SHA256SUMS'), 'w').write('\n'.join(sorted(lines, key=lambda s: s[66:])) + '\n')
