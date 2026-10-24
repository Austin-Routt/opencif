#!/bin/bash
# rc7: trace every shipped mask that a paper's micrograph shows, with the rc7 trace code (the window rule), from the
# rebuilt runs/ of this tree.  Writes provenance/rc7_traces/results/{dinh,gifford,xia}/ and a log per run.
# RC8: the A7 fields (coverage, censored nulls) in every record and caption; the Fig. 4 reproduction is the nested
# reading; rc7's reading B is the variant gifford2014_fig4_pathorder; V13, the taper-direction sign test
# (taper_sign.py), on every Fig. 4 reading: the neighbouring leg at a +-30 um wall window (it must exceed the nested
# reading's ~18-19.7 um wall offsets; RC8 code review, m3), leg 4's array part (control) at +-12 um.
# usage: provenance/rc7_traces/code/run_all_rc7.sh RASTERS [--noscan]
#   RASTERS: the folder with dinh5-000.png, gif5-001.png, gif7-000.png, xia5-000.png (pdfimages -png; not in the repo)
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
R=$(cd "$HERE/../../.." && pwd)
RAS=$1; shift
SCAN=${1:-}
PY="${PYTHON:-python3} -I"
OUT=$R/provenance/rc7_traces/results
mkdir -p $OUT/dinh $OUT/gifford $OUT/xia $OUT/logs
run() {   # name, then the command; skips a run whose log ends with EXIT 0
  local n=$1; shift
  local lg=$OUT/logs/$n.txt
  if grep -q '^EXIT 0' $lg 2>/dev/null; then echo "skip $n (done)"; return; fi
  echo "$(date '+%H:%M:%S') $n"
  "$@" > $lg 2>&1; echo "EXIT $?" >> $lg
  tail -1 $lg
}
mk() { echo "$R/presets/$1.yaml $R/runs/$1/$1.dxf $R/runs/$1/$1.json"; }
# Dinh Fig. 2B: the printed preset and the figure-aligned two-cutoff reading (claims at the same limits: the panels show
# local features both share)
run dinh_rc7 $PY $HERE/run_dinh.py $R $RAS $OUT/dinh rc7 $(mk dinh2024_v1_depicted) "the shipped mask dinh2024_v1_depicted" $SCAN
run dinh_fa  $PY $HERE/run_dinh.py $R $RAS $OUT/dinh fa  $(mk dinh2024_v1_figure_aligned) "the figure-aligned mask dinh2024_v1_figure_aligned" $SCAN
# Gifford Fig. 4 inset: the design preset (reference) and the aligned reproduction (RC8: the nested reading) (claims);
# the rejected readings B (path order) and C (the 51 um start), variants: built here, traced and shown, not anchored
# (their V13 results are)
run gif4_rc7     $PY $HERE/run_gifford.py $R $RAS $OUT/gifford rc7 $(mk gifford2014_fig4_design) "the shipped mask gifford2014_fig4_design" $SCAN
run gif4_aligned $PY $HERE/run_gifford.py $R $RAS $OUT/gifford aligned $(mk gifford2014_fig4_aligned) "the aligned mask gifford2014_fig4_aligned (nested legs, RC8)" $SCAN
V=$R/provenance/rc7_traces/variants
run gif4_pathorder_build $PY $HERE/build_variant.py $R $V/gifford2014_fig4_pathorder.yaml $V/runs/gifford2014_fig4_pathorder
run gif4_pathorder $PY $HERE/run_gifford.py $R $RAS $OUT/gifford pathorder $V/gifford2014_fig4_pathorder.yaml \
    $V/runs/gifford2014_fig4_pathorder/gifford2014_fig4_pathorder.dxf $V/runs/gifford2014_fig4_pathorder/gifford2014_fig4_pathorder.json \
    "the rejected variant gifford2014_fig4_pathorder (rc7's reading B, legs in path order)" $SCAN
run gif4_start51_build $PY $HERE/build_variant.py $R $V/gifford2014_fig4_start51.yaml $V/runs/gifford2014_fig4_start51
run gif4_start51 $PY $HERE/run_gifford.py $R $RAS $OUT/gifford start51 $V/gifford2014_fig4_start51.yaml \
    $V/runs/gifford2014_fig4_start51/gifford2014_fig4_start51.dxf $V/runs/gifford2014_fig4_start51/gifford2014_fig4_start51.json \
    "the rejected variant gifford2014_fig4_start51 (51 um start, nested)" $SCAN
# the 20 um folded-head variant (rc7 final, register INV-2 and BUILD-2): built here, traced and shown, not anchored
run gif4_fold20_build $PY $HERE/build_variant.py $R $V/gifford2014_fig4_fold20.yaml $V/runs/gifford2014_fig4_fold20
run gif4_fold20 $PY $HERE/run_gifford.py $R $RAS $OUT/gifford fold20 $V/gifford2014_fig4_fold20.yaml \
    $V/runs/gifford2014_fig4_fold20/gifford2014_fig4_fold20.dxf $V/runs/gifford2014_fig4_fold20/gifford2014_fig4_fold20.json \
    "the folded-head variant gifford2014_fig4_fold20 (20 um floor, legs in path order)" $SCAN
# V13 (RC8): the taper-direction sign test on every Fig. 4 reading's neighbouring leg (+-30 um) and leg 4 (+-12 um, its
# array part), from each trace's map
v13() {   # tag, the DXF
  run v13_$1_neighbour $PY $HERE/taper_sign.py $RAS/gif5-001.png $OUT/gifford/gifford_results_$1.json $2 neighbour $1_neighbour $OUT/gifford/v13_$1_neighbour.json 30
  run v13_$1_leg4 $PY $HERE/taper_sign.py $RAS/gif5-001.png $OUT/gifford/gifford_results_$1.json $2 leg4 $1_leg4 $OUT/gifford/v13_$1_leg4.json 12 5 3 @array
}
v13 rc7 $R/runs/gifford2014_fig4_design/gifford2014_fig4_design.dxf
v13 aligned $R/runs/gifford2014_fig4_aligned/gifford2014_fig4_aligned.dxf
for vt in pathorder start51 fold20; do v13 $vt $V/runs/gifford2014_fig4_$vt/gifford2014_fig4_$vt.dxf; done
# Gifford Fig. 6: each leg on its own map, and the four legs with their turns on one map, for design and aligned
for t in design aligned; do
  run gif6_${t}        $PY $HERE/run_gifford6.py $R $RAS $OUT/gifford $t $(mk gifford2014_fig6_$t) "the mask gifford2014_fig6_$t" $SCAN
  run gif6_${t}_onemap $PY $HERE/run_gifford6.py $R $RAS $OUT/gifford $t $(mk gifford2014_fig6_$t) "the mask gifford2014_fig6_$t" --one-map --turns $SCAN
done
# Xia Fig. 3: panel g, panels c-f (c and d placed by their own features, e and f at Fig. 3a's box centres), Fig. 3a's
# whole device, for the depicted and aligned presets
for p in "rc7 xia2016_rs7_depicted" "pp xia2016_rs7_aligned"; do
  set -- $p; t=$1; n=$2
  run xia_${t}_fig3a $PY $HERE/xia_fig3a_device.py $RAS/xia5-000.png $R/runs/$n/$n.dxf $R/runs/$n/$n.json $OUT/xia/xia_fig3a_device_$t.json $OUT/xia/xia_fig3a_device_$t.png
  run xia_${t}_g     $PY $HERE/run_xia.py $R $RAS $OUT/xia $t rc6 $(mk $n) "the mask $n" $SCAN
  run xia_${t}_cdef  $PY $HERE/run_xia_more.py $R $RAS $OUT/xia $t cdef $(mk $n) "the mask $n" --fig3a $OUT/xia/xia_fig3a_device_$t.json $SCAN
  run xia_${t}_v12   $PY $HERE/xia_fig3a_device.py $RAS/xia5-000.png $R/runs/$n/$n.dxf $R/runs/$n/$n.json $OUT/xia/xia_fig3a_device_$t.json $OUT/xia/xia_fig3a_device_$t.png \
      $OUT/xia/xia_results_fig3c_$t.json:-:c $OUT/xia/xia_results_fig3d_$t.json:-:d $OUT/xia/xia_results_fig3e_$t.json:-:e $OUT/xia/xia_results_fig3f_$t.json:-:f $OUT/xia/xia_results_$t.json:-:g
done
# Xia panels f and g: the side channels' outer walls against each mask's, on each panel's own registration (register
# INV-4, open)
run xia_side_walls $PY $HERE/xia_side_walls.py $RAS $OUT/xia/xia_side_walls.json \
    f_rc7=$OUT/xia/xia_results_fig3f_rc7.json:$R/runs/xia2016_rs7_depicted/xia2016_rs7_depicted.dxf \
    g_rc7=$OUT/xia/xia_results_rc7.json:$R/runs/xia2016_rs7_depicted/xia2016_rs7_depicted.dxf \
    f_pp=$OUT/xia/xia_results_fig3f_pp.json:$R/runs/xia2016_rs7_aligned/xia2016_rs7_aligned.dxf \
    g_pp=$OUT/xia/xia_results_pp.json:$R/runs/xia2016_rs7_aligned/xia2016_rs7_aligned.dxf
echo "ALL TRACES $(date '+%F %T')"
