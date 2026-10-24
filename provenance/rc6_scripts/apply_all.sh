#!/bin/bash
# rc6 from the rc5 post-review tree (opencif_evaluation.zip of 2026-10-12, code 56e1980d8a064b35): every hashed change,
# in order.  Usage: bash provenance/rc6_scripts/apply_all.sh RASTERS_DIR   (run from the repository root; Python 3.11
# with requirements-tested.txt).  The rasters are the publishers' figures extracted with pdfimages (README).
set -e
RAS=$1
PY=${PYTHON:-python3}
H=provenance/rc6_scripts
mkdir -p opencif/data/r6
$PY -I $H/measure_rc6.py "$PWD" "$RAS" $H/rc6_measurements.json opencif/data/r6/rc6_values.json
$PY $H/apply_rc6_code.py .
$PY $H/apply_rc5_deferred.py .
$PY $H/apply_rc6_e5_guard.py .
$PY $H/apply_rc6_devspecs.py .
$PY $H/apply_rc6_review_fixes.py .                 # the independent review's three code items (R5, R11, R12)
$PY $H/make_figure_anchors_r6.py .
(cd opencif/engine && $PY depicted.py > ../../$H/out/depicted_rc6.out)
$PY -c "from opencif import pipeline; print('code', pipeline.CODE)"
$PY $H/apply_rc6_presets.py .
# the Dinh length computations at the presets' own pitch (about 6 min) and the presets' notes from them
$PY $H/dinh_length_bound.py . presets/dinh2024_v1_depicted.yaml $H/dinh_length_bound_depicted.json > $H/out/dinh_length_bound_depicted.out
$PY $H/dinh_length_bound.py . presets/dinh2024_v1_specified.yaml $H/dinh_length_bound_specified.json > $H/out/dinh_length_bound_specified.out
$PY $H/dinh_taper_alternatives.py . presets/dinh2024_v1_depicted.yaml $H/dinh_taper_alternatives.json $H/dinh_length_bound_depicted.json $H/dinh_length_bound_specified.json > $H/out/dinh_taper_alternatives.out
$PY $H/dinh_schedule_at_F.py . $H/dinh_schedule_at_F.json > $H/out/dinh_schedule_at_F.out
$PY $H/apply_dinh_notes.py .
# Gifford Fig. 6's micrograph (about 3 min) and its two presets
$PY -I $H/measure_gifford_fig6.py "$RAS/gif7-000.png" $H/gifford_fig6_measurement.json > $H/out/measure_gifford_fig6.out
$PY $H/apply_rc6_fig6_presets.py .
# rc6 (after the 2026-10-17 review): the Gifford presets take the layout their figures show (read by the traces on the
# figures alone: provenance/rc6_traces/measurements, results)
$PY $H/measure_gifford_layout.py .
$PY $H/apply_rc6_layout_presets.py .
