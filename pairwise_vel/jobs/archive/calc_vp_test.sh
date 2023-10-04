#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --account=def-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=calc-vp-abacus-planck-z0.70-v-nl-mass-script2-maxs-10
#SBATCH --output=../output/job_logs/abacus/%x-%j.out


echo $INPUT $OUTPUT $FILL_TYPE| ../code/02-calc_vp_halo_llist_pbc_final-var.exe
# echo $INPUT $OUTPUT $FILL_TYPE $NBIN| ../code/02-calc_vp_halo_llist_pbc_final-var.exe
