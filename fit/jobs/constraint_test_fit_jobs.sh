# for run in std1_gl std2_gl std3_gl std1_s8 std2_s8 std1_split std2_split std3_split
for run in std3_gl
  do
    for direction in low high
    do
#    run_name="${run}_${direction}"
    run_name="${run}_${direction}_ctp"
    export INPUT_FILE="/home/mj3chapm/P2/fit/input/fit_inputs/smoothed_emulator/constraint_test/${run_name}.yaml"
    sed "s/RUN_NAME/${run_name}/" mcmc_fit_template.sh > mcmc_fit_run.sh
    sbatch mcmc_fit_run.sh
    sleep 1
    done
  done

# run_name="base"
run_name="base_ctp"
export INPUT_FILE="/home/mj3chapm/P2/fit/input/fit_inputs/smoothed_emulator/constraint_test/${run_name}.yaml"
sed "s/RUN_NAME/${run_name}/" mcmc_fit_template.sh > mcmc_fit_run.sh
sbatch mcmc_fit_run.sh

run_name="std1_s8_low_no_tp"
export INPUT_FILE="/home/mj3chapm/P2/fit/input/fit_inputs/smoothed_emulator/constraint_test/${run_name}.yaml"
sed "s/RUN_NAME/${run_name}/" mcmc_fit_template.sh > mcmc_fit_run.sh
# sbatch mcmc_fit_run.sh
