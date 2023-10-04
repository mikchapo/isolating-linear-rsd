#!/bin/bash
#SBATCH --time=RUN_TIME
#SBATCH --cpus-per-task=NCPUS
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=JOB_NAME
#SBATCH --output=../output/job_logs/OUTPUT_DIR/%x-%j.out


echo $INPUT $OUTPUT $FILL_TYPE| ../code/02-calc_vp_halo_llist_pbc_final-var.exe
