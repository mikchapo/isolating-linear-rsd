#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=analyze-chains-INPUT_NAME
#SBATCH --output=../output/job_logs/analyze_chains/smoothed_emulator/test_clustering/%x-%j.out

module load scipy-stack
source ~/P2/fit/emEnv/bin/activate

python ../code/analyze_chains.py $INPUT_FILE

deactivate
