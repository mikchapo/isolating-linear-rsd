#!/bin/bash
#SBATCH --time=5-0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=125G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=mcmc-fit-RUN_NAME
#SBATCH --output=../output/job_logs/mcmc_fits/%x-%j.out

module load scipy-stack
module load mpi4py
source ~/P2/fit/emEnv/bin/activate

cd ~/P2/fit/code
srun python mcmc_fit.py $INPUT_FILE

deactivate
