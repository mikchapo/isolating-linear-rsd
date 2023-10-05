#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=DD-ang
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

echo $DATA_INPUT $BW_INPUT $OUTPUT_ROOT $PAR_OUTPUT $FIB_OUTPUT $NREG| ../code/01-comp_DD_ang.exe
