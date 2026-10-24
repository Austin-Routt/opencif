"""The whole-device figures (full_device/), end to end and for any tree: writes the configuration
(full_device_config.py), draws each device in its own process (full_device.py; one process per device keeps the memory
under ~3.3 GB, the Xia devices at strip_max_mm 20 and pdf_dpi 300), merges the devices' indexes into
OUT/full_device_index.json and assembles the four shipped devices at one scale (composite_from_strips.py, 4 µm/px).
This is the sequence rc6's full_device/ was made with.

usage: python3 -I run_full_device.py OPENCIF_ROOT RASTERS OUT_DIR [VARIANTS_DIR] [--devices NAME,NAME,...] [--no-composite]
  RASTERS       the publishers' page rasters (README: pdfimages); VARIANTS_DIR the two Xia variant builds
                (build_variant.py; optional: without it the variants are left out)
  OUT_DIR       written: <name>_full.pdf, <name>_strip<i>.png, <name>_zoom_<letter>.png, full_device_composite.png,
                full_device_index.json, config*.json (the configurations, with this machine's paths)"""
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SHIPPED = ['dinh2024_v1_depicted', 'gifford2014_fig4_design', 'gifford2014_fig6_design', 'xia2016_rs7_depicted']


def main(argv):
    flags = [a for a in argv if a.startswith('--')]
    args = [a for a in argv if not a.startswith('--')]
    names_req = None
    if '--devices' in argv:
        names_req = argv[argv.index('--devices') + 1].split(',')
        args = [a for a in args if a != argv[argv.index('--devices') + 1]]
    if len(args) < 3:
        print(__doc__)
        return 1
    root, ras, out = (os.path.abspath(a) for a in args[:3])
    var = os.path.abspath(args[3]) if len(args) > 3 else None
    os.makedirs(out, exist_ok=True)
    traces = os.path.join(root, 'provenance', 'rc6_traces')
    cfg_path = os.path.join(out, 'config.json')
    subprocess.run([sys.executable, '-I', os.path.join(HERE, 'full_device_config.py'), root, ras,
                    os.path.join(traces, 'results'), cfg_path] + ([var] if var else []), check=True)
    cfg = json.load(open(cfg_path))
    names = [d['name'] for d in cfg['devices'] if names_req is None or d['name'] in names_req]
    work = os.path.join(out, '_each')
    merged = dict(config='one process per device (full_device_config.py, then full_device.py on each device\'s own '
                         'config_<name>.json; the Xia devices at strip_max_mm 20 and pdf_dpi 300)', devices={}, failures=[], ok=True)
    for d in cfg['devices']:
        if d['name'] not in names:
            continue
        e = {k: v for k, v in cfg.items() if k != 'devices'}
        e.update(composite=False, devices=[d])
        if d['name'].startswith('xia'):
            e.update(strip_max_mm=20, pdf_dpi=300)
        c1 = os.path.join(out, 'config_%s.json' % d['name'])
        json.dump(e, open(c1, 'w'), indent=1)
        w1 = os.path.join(work, d['name'])
        os.makedirs(w1, exist_ok=True)
        r = subprocess.run([sys.executable, '-I', os.path.join(HERE, 'full_device.py'), c1, w1], cwd=traces)
        idx = json.load(open(os.path.join(w1, 'full_device_index.json')))
        print('%s: exit %d, %d failures' % (d['name'], r.returncode, len(idx.get('failures', []))), flush=True)
        for f in idx['devices'][d['name']]['files']:              # only the files this run's index lists
            shutil.copy(os.path.join(w1, f), out)
        merged['devices'][d['name']] = idx['devices'][d['name']]
        merged['failures'] += idx.get('failures', [])
        merged['ok'] = merged['ok'] and bool(idx.get('ok')) and r.returncode == 0
    json.dump(merged, open(os.path.join(out, 'full_device_index.json'), 'w'), indent=1)
    shutil.rmtree(work)
    if '--no-composite' not in flags and all(n in merged['devices'] for n in SHIPPED):
        subprocess.run([sys.executable, '-I', os.path.join(HERE, 'composite_from_strips.py'), out, '2'] + SHIPPED, check=True)
    print('full_device: %s (%d devices, %d failures)' % ('ok' if merged['ok'] else 'FAILED', len(merged['devices']), len(merged['failures'])))
    return 0 if merged['ok'] else 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
