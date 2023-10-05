#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=mocks-cfs-smooth-test
#SBATCH --output=../output/job_logs/%x-%j.out

source /home/mj3chapm/P2/emulator_model/em_model_env/bin/activate

cd ../code/
python mock_cfs.py AbacusCosmos_1100box 00 rockstar split gamma_l False True

deactivate
