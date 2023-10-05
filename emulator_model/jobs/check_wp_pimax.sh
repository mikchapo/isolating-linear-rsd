#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=check-wp-pimax
#SBATCH --output=../output/job_logs/%x-%j.out

source /home/mj3chapm/P2/emulator_model/em_model_env/bin/activate

cd /home/mj3chapm/P2/emulator_model/code/
python check_wp_pimax.py

deactivate
