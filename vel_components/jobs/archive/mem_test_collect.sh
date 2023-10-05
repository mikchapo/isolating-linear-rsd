#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=mem-test-collect
#SBATCH --output=../output/job_logs/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python mem_test_collect.py
