#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --account=def-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=calc-lin-vel-AbacusCosmos_1100box_planck-00-0-rockstar-True
#SBATCH --output=../output/job_logs/AbacusCosmos_1100box_planck_rockstar/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python calc_lin_vel.py $boxid $halotype $simname $subhalos
