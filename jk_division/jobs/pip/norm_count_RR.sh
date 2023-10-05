#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=norm-count-RR-pip
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

echo $RAND_INPUT $OUTPUT_ROOT $RR_NORM_OUTPUT $NREG| ../code/13-norm_count_RR_pip.exe
