#!/bin/bash
# Regenerate an opencif 1.0.0 working tree from the frozen RC8 extraction:
#   RC8_TREE=<the frozen RC8 tree> PYTHON=<python with the tested versions> bash make_tree.sh OUT_TREE
# The scripts are the ones beside this file (tools_v1/ in the author's work, validation/restructure/ in the release).
set -e
PY=${PYTHON:-python3}
T=$(cd "$(dirname "$0")" && pwd)
RC8=${RC8_TREE:?set RC8_TREE to the frozen RC8 tree}
OUT=$1
$PY -I -B $T/refactor_v1.py $RC8 $OUT $T/kernels.py | head -1
$PY -I -B $T/patch_v1.py $OUT
$PY -I -B $T/patch_audit.py $OUT $T/auditor.py
$PY -I -B $T/patch_cli.py $OUT $T
$PY -I -B $T/patch_names.py $OUT
$PY -I -B $T/patch_lint.py $OUT
$PY -I -B $T/patch_voice.py $OUT $T/voice
$PY -I -B $T/patch_tests.py $OUT
$PY -I -B $T/patch_voice.py $OUT $T/voice tests
$PY -I -B $T/make_glossary_values.py $OUT
$PY -I -B $T/patch_docs.py $OUT $T
$PY -I -B $T/patch_devices.py $OUT $T
$PY -I -B $T/patch_voice.py $OUT $T/voice devices tools
if [ -f $T/voice_docs/MANIFEST.json ]; then $PY -I -B $T/voice_docs.py $OUT $T/voice_docs; fi
# the device kits' source notes and the two generators of documents, whose rewrite changes the text of strings
$PY -I -B $T/voice_docs.py $OUT $T/voice_strings
# the changes made in the tree after the voice pass (patch_final.py names them)
$PY -I -B $T/patch_final.py $OUT $T/final
# the generated documents that the suites require to be current, and the sketch of Figure 1
(cd $OUT && $PY -I -B tools/options_sheet.py)
(cd $OUT && $PY -I -B tools/tutorial_figures.py docs/figures)
