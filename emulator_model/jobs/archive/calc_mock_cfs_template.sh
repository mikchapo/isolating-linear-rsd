#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=calc-mocks-cfs-SIMNAME-BOXID-HALOTYPE-SCALING_METHOD-SCALE_VARIABLE-SMOOTHED_VEL-N_GRID_SMOOTH_TYPE_R_SMOOTH
#SBATCH --output=../output/job_logs/%x-%j.out

source /home/mj3chapm/P2/emulator_model/em_model_env/bin/activate

cd /home/mj3chapm/P2/emulator_model/code/
python calc_mock_cfs.py $simname $boxid $halotype $scaling_method $scale_variable $load_dvs $smoothed_vel $N_grid $smooth_type $R_smooth

deactivate
