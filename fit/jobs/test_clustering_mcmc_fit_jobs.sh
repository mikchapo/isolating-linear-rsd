# for ID in {0..9}
for ID in 6
  do
#  for VOL_FACTOR in eBOSS 5 20 60
  for VOL_FACTOR in eBOSS 5
    do
    for RUN in default rm_small rm_large large_only intermediate_only
      do
      RUN_NAME="test_hod_${ID}_vol-factor-${VOL_FACTOR}_${RUN}"
      export INPUT_FILE="/home/mj3chapm/P2/fit/input/fit_inputs/smoothed_emulator/test_hod/${RUN_NAME}.yaml"
      sed "s/RUN_NAME/${RUN_NAME}/" test_clustering_mcmc_fit_template.sh > test_clustering_mcmc_fit_run.sh
      sbatch test_clustering_mcmc_fit_run.sh
      sleep 1
      done
    done
  done
