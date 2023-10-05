#!/bin/bash
#SBATCH --time=01:20:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=mocks-cfs-slv-planck-00-0-rockstar-halos-split-gamma_l
#SBATCH --output=../output/job_logs/%x-%j.out

source /home/mj3chapm/P2/emulator_model/em_model_env/bin/activate

python ../code/mock_cfs_slv.py rockstar split gamma_l 1. 1. False

deactivate
