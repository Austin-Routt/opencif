"""Assemble a device package from this repository and run it: the photomasks of a published device built with this
release, registered on the paper's micrographs, measured and rendered (devices/README.md).

usage: python3 tools/device_package.py DEVICE --rasters DIR --out DIR [STEP ...]
  DEVICE    dinh2024_v1 | gifford2014_fig4 | gifford2014_fig6 | xia2016_rs7
  --rasters the folder holding the publishers' figures the device needs (devices/<DEVICE>/MANIFEST.json names each and
            its sha256; they are not in this repository)
  --out     where the package is assembled: OUT/opencif_<DEVICE>/
  STEP      check, build, trace, render ... (default: every step of the device's script)

The package is a clean room: the repository's files the device's MANIFEST.json lists are copied into
OUT/opencif_<DEVICE>/tree/ (each checked against its sha256), the kit's own files beside them, the figures from
--rasters (checked by sha256), and the device's script runs there in its own process.  The script checks the code
stamp, every input, each photomask's sha256 and every trace record against the release's records, and stops with
exit 2 at the first departure."""
import argparse
import hashlib
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def assemble(name, rasters, out):
    man = json.load(open(os.path.join(ROOT, 'devices', name, 'MANIFEST.json')))
    pkg = os.path.join(os.path.abspath(out), 'opencif_' + name)
    for rel, src in man['files'].items():
        if rel.startswith('/') or '\\' in rel or ':' in rel or '..' in rel.split('/'):
            raise SystemExit('unsafe path in the %s kit: %r' % (name, rel))
        target = os.path.join(pkg, *rel.split('/'))
        os.makedirs(os.path.dirname(target), exist_ok=True)
        if 'raster' in src:
            cand = os.path.join(rasters, os.path.basename(rel)) if rasters else None
            if not cand or not os.path.exists(cand):
                raise SystemExit('%s needs the publisher\'s figure %s (sha256 %s...) in --rasters'
                                 % (name, os.path.basename(rel), src['raster'][:16]))
            data = open(cand, 'rb').read()
            if sha256_bytes(data) != src['raster']:
                raise SystemExit('%s: %s is not the expected image (sha256 %s, expected %s)'
                                 % (name, cand, sha256_bytes(data)[:16], src['raster'][:16]))
        else:
            data = open(os.path.join(ROOT, *src['repo'].split('/')), 'rb').read()
            if sha256_bytes(data) != src['sha256']:
                raise SystemExit('%s: %s (from %s) is not the file the kit lists' % (name, rel, src['repo']))
        if not (os.path.exists(target) and open(target, 'rb').read() == data):
            with open(target, 'wb') as fh:
                fh.write(data)
    # The package keeps its MANIFEST.json, against which its check step verifies every file again.
    with open(os.path.join(pkg, 'MANIFEST.json'), 'w') as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    return pkg, man['script']


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('device')
    ap.add_argument('--rasters', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('steps', nargs='*')
    a = ap.parse_args(argv)
    pkg, script = assemble(a.device, a.rasters, a.out)
    print('assembled %s' % pkg, flush=True)
    return subprocess.call([sys.executable, '-I', '-B', script] + list(a.steps), cwd=pkg)


if __name__ == '__main__':
    sys.exit(main())
