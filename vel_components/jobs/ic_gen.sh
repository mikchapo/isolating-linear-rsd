#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --account=def-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=ic-gen-AbacusCosmos_1100box_planck-00-19
#SBATCH --output=../output/job_logs/%x-%j.out

cd ../zeldovich-PLT/
./zeldovich /home/mj3chapm/scratch/abacus/${simname}_products/${simname}_${boxid}_products/${simname}_${boxid}_rockstar_halos/info/abacus_ngc_nplt.par
