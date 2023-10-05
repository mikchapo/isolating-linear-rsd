#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=DR-mps-pip-ang
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $BW_INPUT $INPUT_ROOT $DR_AUW_INPUT $OUTPUT_ROOT $DR_OUTPUT $NREG| ../code/06-comp_DR_mps_pip.exe
