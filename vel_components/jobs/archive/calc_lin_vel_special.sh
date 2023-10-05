#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --account=def-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=calc-lin-vel-AbacusCosmos_1100box-planck-rockstar
#SBATCH --output=../output/job_logs/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python calc_lin_vel.py planck rockstar AbacusCosmos_1100box
