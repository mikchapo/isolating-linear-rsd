#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --mem=125G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=check-offset
#SBATCH --output=../output/job_logs/%x-%j.out

module load python mpi4py
source ../vel_comp_env/bin/activate
cd ../code/
boxid=00-0
halotype=rockstar
simname=AbacusCosmos_1100box_planck
mpirun -np 12 python check_offset.py $boxid $halotype $simname
