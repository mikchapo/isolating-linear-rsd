#!/bin/bash
#SBATCH --time=RUN_TIME
#SBATCH --cpus-per-task=NCPUS
#SBATCH --mem-per-cpu=2G
#SBATCH --account=def-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=JOB_NAME
#SBATCH --output=../output/job_logs/OUTPUT_DIR/%x-%j.out


echo $INPUT $OUTPUT $FILL_TYPE| ../code/SCRIPT
# echo $INPUT $OUTPUT $FILL_TYPE $NBIN| ../code/SCRIPT
