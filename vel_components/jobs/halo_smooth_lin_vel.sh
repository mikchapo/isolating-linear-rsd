#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --begin=now+10hour
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=140G
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=END
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=halo-smooth-lin-vel-AbacusCosmos_1100box_planck-00-19-rockstar-1100-tophat-5-v3
#SBATCH --output=../output/job_logs/vel_smoothing/v3/test/%x-%j.out

source ../vel_comp_env/bin/activate
cd ../code/
python halo_smooth_lin_vel.py $simname $boxid $halo_type $N_grid $smooth_type $R_smooth $subhalos
