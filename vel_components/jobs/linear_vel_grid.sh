#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=END
#SBATCH --job-name=linear-vel-grid-AbacusCosmos_1100box_planck-00-19-1375-v3
#SBATCH --output=../output/job_logs/vel_smoothing/v3/test/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python linear_vel_grid.py $boxid $simname $N_grid
