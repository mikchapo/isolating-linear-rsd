#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=norm-count-DR-pip
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $BW_INPUT $OUTPUT_ROOT $DR_NORM_OUTPUT $NREG| ../code/12-norm_count_DR_pip.exe
