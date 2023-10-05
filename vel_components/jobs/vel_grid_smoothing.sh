#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --begin=now+5hour
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=80G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=END
#SBATCH --job-name=vel-grid-smoothing-AbacusCosmos_1100box_planck-00-19-1100-gaussian-2-False-v3
#SBATCH --output=../output/job_logs/vel_smoothing/v3/test/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python vel_grid_smoothing.py $boxid $simname $N_grid $smooth_type $R_smooth $fft_smooth
