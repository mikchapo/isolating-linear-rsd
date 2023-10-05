#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=norm-count-DR-ang
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $BW_INPUT $OUTPUT_ROOT $NORM_OUTPUT $NREG| ../code/04-norm_count_DR_ang.exe
