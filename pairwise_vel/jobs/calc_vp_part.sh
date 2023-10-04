#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=calc-vp-part-z49-ic-1000000-ngc-v6
#SBATCH --output=../output/job_logs/%x-%j.out

cd ../code/
INPUT=/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_sample_z49.0_1000000_ic_no_growth_corr_v6.dat
OUTPUT=/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z49.0_1000000_ic_no_growth_corr_vp_v6.dat
echo $INPUT $OUTPUT| ./04-calc_vp_part_llist_pbc_var.exe
