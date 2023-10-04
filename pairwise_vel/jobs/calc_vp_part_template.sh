#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=calc-vp-part-abacus-training-BOXID
#SBATCH --output=../output/job_logs/abacus/training_part_test/%x-%j.out

cd ../code/
echo $INPUT $OUTPUT| ./04-calc_vp_part_llist_pbc_var.exe
