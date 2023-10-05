#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=norm-count-DD-ang
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

echo $DATA_INPUT $BW_INPUT $OUTPUT_ROOT $NORM_OUTPUT $NREG| ../code/03-norm_count_DD_ang.exe
