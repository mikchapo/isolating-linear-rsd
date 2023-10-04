for ID in {0..9}
  do
  for RUN in vol-factor-20_default vol-factor-20_scales vol-factors
    do
    INPUT_NAME="test_hod_${ID}_${RUN}"
    export INPUT_FILE="/home/mj3chapm/P2/fit/input/tp_inputs/smoothed_emulator/test_hod/${INPUT_NAME}.yaml"
    sed "s/INPUT_NAME/${INPUT_NAME}/" test_clustering_analyze_chains_template.sh > test_clustering_analyze_chains_run.sh
#    sbatch test_clustering_analyze_chains_run.sh
#    sleep 1
    done
  done

for VOL_FACTOR in eBOSS 5 60
  do
  for RUN in default scales
    do
    INPUT_NAME="test_hod_0_vol-factor-${VOL_FACTOR}_${RUN}"
    export INPUT_FILE="/home/mj3chapm/P2/fit/input/tp_inputs/smoothed_emulator/test_hod/${INPUT_NAME}.yaml"
    sed "s/INPUT_NAME/${INPUT_NAME}/" test_clustering_analyze_chains_template.sh > test_clustering_analyze_chains_run.sh
    sbatch test_clustering_analyze_chains_run.sh
    sleep 1
    done
  done
