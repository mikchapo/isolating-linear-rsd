#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=calc-vp-test-00-19-v-tot-mass
#SBATCH --output=../output/job_logs/abacus/%x-%j.out


echo $INPUT $OUTPUT $FILL_TYPE| ../code/02-calc_vp_halo_llist_pbc_final-var.exe
