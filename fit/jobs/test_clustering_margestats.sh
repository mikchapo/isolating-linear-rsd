#!/bin/bash
#SBATCH --array=0-9
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=test-clustering-margestats
#SBATCH --output=../output/job_logs/analyze_chains/smoothed_emulator/test_clustering/%x-%j.out

RUNNAME=test_hod_${SLURM_ARRAY_TASK_ID}_vol-factor-20_default

cd ~/projects/rrg-wperciva/mj3chapm/P2/chains/smoothed_emulator/test_clustering/${RUNNAME}
getdist --ignore_rows=40000 ./$RUNNAME

