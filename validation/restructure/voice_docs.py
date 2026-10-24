"""voice_docs.py: put the documents of opencif 1.0.0 that RC8 wrote in the author's voice.

Usage: python -I voice_docs.py TREE VOICE_DOCS_DIR

VOICE_DOCS_DIR holds the rewritten documents under their paths in the tree, and MANIFEST.json, which gives for each
the sha256 of the document it was rewritten from (the document as patch_docs.py leaves it). A rewritten document
replaces the tree's only when the tree's document is still that one, so a later edit of patch_docs.py cannot be lost
under an older rewrite. Applied after patch_docs.py.
"""
import hashlib
import json
import os
import shutil
import sys

tree, vdir = sys.argv[1], sys.argv[2]
manifest = json.load(open(os.path.join(vdir, 'MANIFEST.json'), encoding='utf-8'))
stale = []
for rel, pre in sorted(manifest.items()):
    dst = os.path.join(tree, *rel.split('/'))
    with open(dst, 'rb') as fh:
        if hashlib.sha256(fh.read()).hexdigest() != pre:
            stale.append(rel)
            continue
    shutil.copy(os.path.join(vdir, *rel.split('/')), dst)
if stale:
    raise SystemExit('voice_docs: rewritten from an older version, so not applied: %s' % stale)
print('voice_docs: %d documents in the author\'s voice' % len(manifest))
