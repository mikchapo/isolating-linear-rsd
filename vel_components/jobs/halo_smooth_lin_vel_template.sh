#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --begin=now+10hour
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=140G
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=END
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=halo-smooth-lin-vel-SIMNAME-BOXID-HALO_TYPE-N_GRID-SMOOTH_TYPE-R_SMOOTH-v3
#SBATCH --output=../output/job_logs/vel_smoothing/v3/OUTPUT_DIR/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python halo_smooth_lin_vel.py $simname $boxid $halo_type $N_grid $smooth_type $R_smooth $subhalos
