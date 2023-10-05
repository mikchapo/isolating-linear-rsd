#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=linear-vel-grid-test
#SBATCH --output=../output/job_logs/vel_smoothing/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python linear_vel_grid_test.py 00 AbacusCosmos_1100box 1100
