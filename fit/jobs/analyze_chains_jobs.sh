# for input_name in eboss_fiducial eboss_equal_gammas eboss_scales eboss_probes
for input_name in eboss_scales
  do
  export INPUT_FILE="/home/mj3chapm/P2/fit/input/tp_inputs/smoothed_emulator/jk-corr/${input_name}.yaml"
  sed "s/INPUT_NAME/${input_name}/" analyze_chains_template.sh > analyze_chains_run.sh
  sbatch analyze_chains_run.sh
  sleep 1
  done
