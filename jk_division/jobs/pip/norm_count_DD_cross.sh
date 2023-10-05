#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=norm-count-DD-cross-pip
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

echo $DATA_INPUT_NGC $DATA_INPUT_SGC $BW_INPUT_NGC $BW_INPUT_SGC $OUTPUT_ROOT $DD_CROSS_NORM_OUTPUT $NREG| ../code/14-norm_count_DD_pip_cross.exe
