#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=calc-ang-weight
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

python ../code/calc_ang_weight.py $INPUT_ROOT $PAR_INPUT $FIB_INPUT $NORM_INPUT $OUTPUT_ROOT $WEIGHT_OUTPUT $NREG
