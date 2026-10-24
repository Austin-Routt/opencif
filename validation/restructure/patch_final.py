"""patch_final.py: install the changes made in the tree after the voice pass (make_tree.sh, before the generated
documents).

Usage: python -I patch_final.py TREE FINAL_DIR

FINAL_DIR holds the final version of each file that was edited after the voice pass, at its place in the tree, and
MANIFEST.json, which gives for each file the sha256 of the version that the steps before this one write, or null for a
file that they do not write. A file is replaced only when the tree's version has that sha256, and added only when it
is absent, so this step changes exactly the files it names (a diff of the tree before and after this step shows each
change). The changes are:
  - tests/test_auditor.py, part 4: the ruler's corners on a synthetic taper that changes slope at one gap
    (CHANGELOG item 12, with the ruler's wall_corners, whose source is tools_v1/ruler.py);
  - tests/test_composition.py: a case that declares Q1 or Q2 may fail the physics audit only on the same questions
    (CHANGELOG item 6);
  - the provenance of rc6's traces and one figure measurement: absolute folders of the author's machine replaced by
    paths relative to the working folder;
  - docs/ and foundation/ERRATA_NEXT_UNFREEZE.md: the changelog's items on the ruler, the m14 record, the audit rule and
    erratum D107, the validation and evaluation pages' tables, and the limitations and audit pages, each still with the
    placeholders that tools_v1/fill_docs.py fills from the final evidence;
  - after the two independent reviews of the 1.0.0 code (CHANGELOG item 14): the documents, README.md, NOTICE.md, the
    workflow (.github/workflows/tests.yml), the options sheet's generator (tools/options_sheet.py), and
    provenance/rc4_scripts/ and rc5_scripts/ with the release candidates' folders written relative to the working
    folder;
  - after the device kits' first run (CHANGELOG item 14): the five trace scripts of provenance/rc7_traces/code/ and
    their copies in rc6_traces/code/, whose one-element tuples of folders read (R,) again;
  - after the independent re-audit of docs/ALIGNMENT.md (CHANGELOG item 14): the matrix (row AA16 relabelled, the new
    row S3, and ten smaller corrections), erratum D108 in foundation/ERRATA_NEXT_UNFREEZE.md, and DISSERTATION_NOTES.
"""
import hashlib
import json
import os
import shutil
import sys

TREE, FINAL = sys.argv[1], sys.argv[2]
man = json.load(open(os.path.join(FINAL, 'MANIFEST.json'), encoding='utf-8'))


def sha256(p):
    with open(p, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


bad = []
for rel, expect in sorted(man.items()):
    dst = os.path.join(TREE, *rel.split('/'))
    if (expect is None and os.path.exists(dst)) or (expect is not None and (not os.path.exists(dst) or
                                                                            sha256(dst) != expect)):
        bad.append(rel)
if bad:
    raise SystemExit('patch_final: these files are not the versions expected before this step: %s' % ', '.join(bad))
for rel in sorted(man):
    dst = os.path.join(TREE, *rel.split('/'))
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(os.path.join(FINAL, *rel.split('/')), dst)
print('patch_final: %d files replaced and %d added, each checked first (by its sha256, or as absent)'
      % (sum(1 for v in man.values() if v is not None), sum(1 for v in man.values() if v is None)))
