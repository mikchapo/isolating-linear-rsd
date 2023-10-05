#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=80G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=test-vel-grid
#SBATCH --output=../output/job_logs/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python test_vel_grid.py
