# for input_name in base std1_gl_low std1_gl_high std2_gl_low std2_gl_high std3_gl_low std3_gl_high std2_om_low std2_s8_low std2_s8_high std2_split_low std2_split_high
for input_name in std2_om_high base_ctp std3_gl_low_ctp std3_gl_high_ctp
  do
  export INPUT_FILE="/home/mj3chapm/P2/fit/input/tp_inputs/smoothed_emulator/constraint_test/${input_name}.yaml"
  sed "s/INPUT_NAME/${input_name}/" analyze_chains_template.sh > analyze_chains_run.sh
  sbatch analyze_chains_run.sh
  sleep 1
  done
