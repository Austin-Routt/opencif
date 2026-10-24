"""make_device_manifests.py: write the MANIFEST.json of each device kit of opencif 1.0.0.

Usage: python -I make_device_manifests.py TREE RC8_MANIFESTS_DIR

A kit's MANIFEST.json names every file of its package with its source and sha256, so that tools/device_package.py can
assemble the package from the repository and the package can check itself. The kits of RC8 carried a copy of RC8's
code under rc8/. The kits of 1.0.0 carry the repository's own package under tree/ instead: every file of opencif/,
requirements-tested.txt, validation/RELEASE.json (the release's code stamp and photomasks), and the same provenance
files and tools that the RC8 kit took, at their places in the tree. The kit's own files come from devices/, with the
shared devicekit.py from devices/_shared/, and the publishers' figures stay as sha256 entries that --rasters must
supply. Each kit's reference/EXPECT.json is then given the release's code stamp, which the device scripts compare with
the package's, and the sha256 of every input that its device script checks, listed under the input's place in the 1.0.0
kit (RC8's rc8/ is tree/) and taken from the files just listed, after which its own sha256 enters the manifest.
"""
import hashlib
import json
import os
import sys

TREE, RC8_DIR = sys.argv[1], sys.argv[2]
sys.path.insert(0, TREE)
release = json.load(open(os.path.join(TREE, 'validation', 'RELEASE.json')))
STAMP = release['code_stamp']


def sha256(rel):
    with open(os.path.join(TREE, *rel.split('/')), 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


for device in sorted(d for d in os.listdir(RC8_DIR) if d.endswith('.json')):
    name = device[:-len('.json')]
    old = json.load(open(os.path.join(RC8_DIR, device)))
    exp_rel = 'devices/%s/reference/EXPECT.json' % name
    files = {}
    for rel, src in old['files'].items():
        if 'raster' in src:
            files[rel] = dict(raster=src['raster'])
        elif rel == 'devicekit.py':
            files[rel] = dict(repo='devices/_shared/devicekit.py')
        elif rel.startswith('rc8/'):
            inner = rel[len('rc8/'):]
            if inner == 'TRACE_CODE_SHA256.txt' or inner.startswith('opencif/'):
                continue                                  # replaced by MANIFEST.json and the whole package below
            files['tree/' + inner] = dict(repo=src.get('tree') or inner)
        elif 'bundle' in src:
            files[rel] = dict(repo=src['bundle'])
        else:
            files[rel] = dict(repo=src['tree'])
    for root, dirs, fns in os.walk(os.path.join(TREE, 'opencif')):
        dirs[:] = sorted(d for d in dirs if d != '__pycache__')
        for fn in sorted(fns):
            if fn.endswith('.pyc'):
                continue
            rel = os.path.relpath(os.path.join(root, fn), TREE).replace(os.sep, '/')
            files['tree/' + rel] = dict(repo=rel)
    for rel in ('requirements-tested.txt', 'validation/RELEASE.json'):
        files['tree/' + rel] = dict(repo=rel)
    for rel, src in files.items():
        if 'repo' in src and src['repo'] != exp_rel:
            src['sha256'] = sha256(src['repo'])
    # The kit's expectations: the release's code stamp, and the inputs that the device script checks, each under its
    # place in this kit with the sha256 of the file that the kit holds there. An input missing from the kit stops here.
    exp = json.load(open(os.path.join(TREE, *exp_rel.split('/'))))
    exp['code_stamp'] = STAMP
    inputs = {}
    for rel in exp['inputs_sha256']:
        new_rel = 'tree/' + rel[len('rc8/'):] if rel.startswith('rc8/') else rel
        if new_rel not in files:
            raise SystemExit('%s: the input %s is not in the kit' % (name, new_rel))
        inputs[new_rel] = files[new_rel].get('sha256') or files[new_rel]['raster']
    exp['inputs_sha256'] = inputs
    with open(os.path.join(TREE, *exp_rel.split('/')), 'w') as fh:
        json.dump(exp, fh, indent=1)
        fh.write('\n')
    files['reference/EXPECT.json']['sha256'] = sha256(exp_rel)
    man = dict(device=name, script=old['script'], code_stamp=STAMP,
               note='assembled by tools/device_package.py: tree/ from this repository, the kit\'s own files from '
                    'devices/, and the figures from --rasters (each checked by its sha256)',
               files=dict(sorted(files.items())))
    out = os.path.join(TREE, 'devices', name, 'MANIFEST.json')
    with open(out, 'w') as fh:
        json.dump(man, fh, indent=1)
        fh.write('\n')
    print('%s: %d files (%d from tree/, %d publishers\' figures)' % (
        name, len(files), sum(1 for r in files if r.startswith('tree/')), sum(1 for s in files.values() if 'raster' in s)))
