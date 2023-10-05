#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=jk-corr-funcs-test
#SBATCH --output=../output/job_logs/%x-%j.out

source ~/P2/jk_division/jk_div_env/bin/activate

cd ~/P2/jk_division/code/
python jk_corr_funcs.py "../output_v1/data_products/jk_pc/" "../output/data_products/eBOSS_pip+ang_jk-corr_200/removed_region_" ".dat" 200

deactivate
