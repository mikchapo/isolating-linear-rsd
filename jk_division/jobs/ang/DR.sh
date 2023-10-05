#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=DR-ang
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

echo $DATA_INPUT $RAND_INPUT $BW_INPUT $OUTPUT_ROOT $PAR_OUTPUT $FIB_OUTPUT $NREG| ../code/02-comp_DR_ang.exe
