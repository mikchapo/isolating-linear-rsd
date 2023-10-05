#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=DR-wp-pip-ang
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $BW_INPUT $INPUT_ROOT $DR_AUW_INPUT $OUTPUT_ROOT $DR_OUTPUT $NREG| ../code/09-comp_DR_wp_pip.exe
