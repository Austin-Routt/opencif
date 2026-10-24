S=scratch; T=$S/rc6patch/tree3; cd $S/rc6patch
mkdir -p final; rm -f final/*; : > final/status.txt
run() { t0=$(date +%s); python3 run_forced.py $T "${@:2}" > final/out_$1.txt 2>&1; echo "$1 exit $? $(( $(date +%s)-t0 )) s" >> final/status.txt; }
M02=examples/matrix/m02_linear_shear_arms_filter_teardrop.yaml
M04=examples/matrix/m04_linear_shear_pack_plain_ring.yaml
M12=examples/matrix/m12_paperclip_shear_arms_plain_ring.yaml
run m02_forced_strict_check 1.5 check $M02
run m04_forced_strict_check 1.5 check $M04
run m12_forced_strict_check 1.5 check $M12
run m02_forced_evaluate_check 1.5 check $S/rc6patch/specs/m02_linear_shear_arms_filter_teardrop_evaluate.yaml
run m12_forced_evaluate_check 1.5 check $S/rc6patch/specs/m12_paperclip_shear_arms_plain_ring_evaluate.yaml
run m02_none_check none check $M02
run m04_none_check none check $M04
run m12_none_check none check $M12
run m02_none_build none build $M02 --out $S/rc6patch/final_build_m02_none
run m02_forced_evaluate_build 1.5 build $S/rc6patch/specs/m02_linear_shear_arms_filter_teardrop_evaluate.yaml --out $S/rc6patch/final_build_m02_forced_eval
echo DONE >> final/status.txt
