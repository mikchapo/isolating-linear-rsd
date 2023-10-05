#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=test-vis-plots
#SBATCH --output=../output/job_logs/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python test_vis_plots.py
