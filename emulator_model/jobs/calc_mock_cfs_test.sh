#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=All
#SBATCH --job-name=calc-mocks-cfs-AbacusCosmos_1100box_planck-00-19-rockstar-split-gamma_n-True-1100_gaussian_2
#SBATCH --output=../output/job_logs/test/%x-%j.out

source /home/mj3chapm/P2/emulator_model/em_model_env/bin/activate

cd /home/mj3chapm/P2/emulator_model/code/
python mock_cfs_test.py $simname $boxid $halotype $scaling_method $scale_variable $smoothed_vel $N_grid $smooth_type $R_smooth $subhalos

deactivate
