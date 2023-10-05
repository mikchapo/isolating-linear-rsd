#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=check-smoothed-vel-dist
#SBATCH --output=../output/job_logs/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python check_smoothed_vel_dist.py
