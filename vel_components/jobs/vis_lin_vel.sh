#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=vis-lin-vel
#SBATCH --output=../output/job_logs/%x-%j.out

module load scipy-stack/2022a
source ../vel_comp_env/bin/activate

cd ../code/
# python vis_lin_vel.py 100. 100. 75.
python vis_lin_vel.py 50. 145. 115.
