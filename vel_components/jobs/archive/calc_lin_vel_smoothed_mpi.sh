#!/bin/bash
#SBATCH --time=5-0
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --mem=125G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=calc-lin-vel-smoothed-mpi-AbacusCosmos_1100box_planck-00-0-rockstar
#SBATCH --output=../output/job_logs/AbacusCosmos_1100box_planck_rockstar/%x-%j.out

module load python mpi4py
source ../vel_comp_env/bin/activate
cd ../code/
boxid=00-0
halotype=rockstar
simname=AbacusCosmos_1100box_planck
mpirun -np 32 python calc_lin_vel_smoothed_mpi.py $boxid $halotype $simname
