#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=def-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=mpi-test
#SBATCH --output=../output/job_logs/%x-%j.out

source ../vel_comp_env/bin/activate
module load python mpi4py
cd ../code/
mpirun -np 3 python3 mpi_test.py
