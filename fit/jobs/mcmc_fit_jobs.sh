# for run_name in uchuu_sham uchuu_sham_rm_small uchuu_sham_rm_large uchuu_sham_intermediate_only uchuu_sham_large_only
# for run_name in eboss_fiducial eboss_rm_small eboss_rm_intermediate eboss_rm_large eboss_large_only eboss_mono+quad eboss_mono+wp eboss_fixed_gamma_n eboss_fixed_gamma_l eboss_equal_gammas uchuu_sham
# for run_name in eboss_fiducial eboss_rm_small eboss_rm_intermediate eboss_rm_large eboss_large_only
for run_name in eboss_equal_gammas_large_only
  do
  export INPUT_FILE="/home/mj3chapm/P2/fit/input/fit_inputs/smoothed_emulator/jk-corr/${run_name}.yaml"
  sed "s/RUN_NAME/${run_name}-jk-corr/" mcmc_fit_template.sh > mcmc_fit_run.sh
  sbatch mcmc_fit_run.sh
  sleep 1
  done
