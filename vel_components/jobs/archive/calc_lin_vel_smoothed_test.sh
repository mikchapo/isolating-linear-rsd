#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=calc-lin-vel-smoothed-test
#SBATCH --output=../output/job_logs/AbacusCosmos_1100box_planck_rockstar/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python calc_lin_vel_smoothed_test.py
