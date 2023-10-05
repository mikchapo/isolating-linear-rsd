#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=8G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=DD-wp-pip-ang
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

echo $DATA_INPUT $BW_INPUT $INPUT_ROOT $DD_AUW_INPUT $OUTPUT_ROOT $DD_OUTPUT $NREG| ../code/08-comp_DD_wp_pip.exe
