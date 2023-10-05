#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=vis-all-part-vel
#SBATCH --output=../output/job_logs/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python vis_lin_vel.py
