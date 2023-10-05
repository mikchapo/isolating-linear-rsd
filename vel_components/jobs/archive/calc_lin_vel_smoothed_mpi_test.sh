#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --ntasks=4
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=calc-lin-vel-smoothed-mpi-AbacusCosmos_1100box_planck-00-0-rockstar-test
#SBATCH --output=../output/job_logs/AbacusCosmos_1100box_planck_rockstar/%x-%j.out

module load python mpi4py
source ../vel_comp_env/bin/activate
cd ../code/
mpirun -np 4 python calc_lin_vel_smoothed_mpi_test.py
