#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=def-wperciva
#SBATCH --job-name=avoid-purge-SIMNAME-BOXID
#SBATCH --output=../output/job_logs/%x-%j.out

cd ~/scratch/abacus/${simname}_products/${simname}_${boxid}_products/
~/scratch/avoid_purge.sh
