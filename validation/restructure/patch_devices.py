"""patch_devices.py: turn RC8's device kits into the kits of opencif 1.0.0.

Usage: python -I patch_devices.py TREE TOOLS_V1

The four kits share one devicekit.py (devices/_shared/devicekit.py), which runs the repository's own package rather
than a copy of RC8's code, and they are assembled by tools/device_package.py. RC8's hash list of the trace code is
replaced by each package's MANIFEST.json, which make_device_manifests.py writes once the release's evidence exists.
devices/README.md describes the kits. Applied after patch_docs.py.
"""
import os
import shutil
import sys

TREE, HERE = sys.argv[1], sys.argv[2]
DEV = os.path.join(TREE, 'devices')
kits = sorted(d for d in os.listdir(DEV) if os.path.isdir(os.path.join(DEV, d)) and not d.startswith('_'))
for k in kits:
    p = os.path.join(DEV, k, 'devicekit.py')
    if os.path.exists(p):
        os.remove(p)
shutil.rmtree(os.path.join(DEV, '_shared', 'rc8'), ignore_errors=True)
shutil.copy(os.path.join(HERE, 'devices_kit', 'devicekit.py'), os.path.join(DEV, '_shared', 'devicekit.py'))
shutil.copy(os.path.join(HERE, 'device_package.py'), os.path.join(TREE, 'tools', 'device_package.py'))
if os.path.exists(os.path.join(HERE, 'devices_kit', 'README.md')):
    shutil.copy(os.path.join(HERE, 'devices_kit', 'README.md'), os.path.join(DEV, 'README.md'))
# no em dashes in what the release writes: the provenance page's title and its empty cell
fp = os.path.join(TREE, 'tools', 'dissertation_figures.py')
s = open(fp, encoding='utf-8').read()
for old, new in (("L = ['# %s \u2014 provenance' % name, '',", "L = ['# %s: provenance' % name, '',"),
                 ("if ((v.get('passage') or {}).get('classes') or {}).get('straight side channel') else '\u2014'),",
                  "if ((v.get('passage') or {}).get('classes') or {}).get('straight side channel') else 'not measured'),")):
    if s.count(old) != 1:
        raise SystemExit('patch_devices: tools/dissertation_figures.py: expected 1, found %d: %r' % (s.count(old), old))
    s = s.replace(old, new)
open(fp, 'w', encoding='utf-8').write(s)
print('patch_devices: %d kits share devices/_shared/devicekit.py; tools/device_package.py assembles them' % len(kits))
