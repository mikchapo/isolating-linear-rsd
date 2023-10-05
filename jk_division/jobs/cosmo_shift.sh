#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=cosmo-shift
#SBATCH --output=../output/%x-%j.out

module load python/3.6

source ~/env/bin/activate

python ../code/cosmo_shift.py

deactivate
