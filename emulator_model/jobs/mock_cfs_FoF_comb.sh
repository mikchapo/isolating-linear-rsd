#!/bin/bash
#SBATCH --time=01:20:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=mocks-cfs-planck-00-0-FoF-halos-combined
#SBATCH --output=../output/job_logs/%x-%j.out

source /home/mj3chapm/RSD/uchuu_mocks/uchuu_mock_env_3.7/bin/activate

python ../code/mock_cfs.py FoF combined gamma_f False

deactivate
