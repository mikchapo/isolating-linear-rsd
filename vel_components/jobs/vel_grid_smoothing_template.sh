#!/bin/bash
#SBATCH --time=HOURS:00:00
#SBATCH --begin=now+5hour
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=MEMORY
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=END
#SBATCH --job-name=vel-grid-smoothing-SIMNAME-BOXID-N_GRID-SMOOTH_TYPE-R_SMOOTH-FFT_SMOOTH-v3
#SBATCH --output=../output/job_logs/vel_smoothing/v3/OUTPUT_DIR/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python vel_grid_smoothing.py $boxid $simname $N_grid $smooth_type $R_smooth $fft_smooth
