#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=32
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=RR-mps-pip
#SBATCH --output=../output/job_logs/jk_pc/%x-%j.out

echo $RAND_INPUT $OUTPUT_ROOT $RR_OUTPUT $NREG| ../code/07-comp_RR_mps.exe
