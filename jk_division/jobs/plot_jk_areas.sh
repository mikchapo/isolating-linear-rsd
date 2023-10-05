#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-wperciva
#SBATCH --job-name=plot-jk-area-data
#SBATCH --output=../output/%x-%j.out

module load python/3.6

source ~/env/bin/activate

python ../code/plot_jk_areas.py ../../corrfunc/multipole_pair_counts/data/eBOSS_LRG_NGC_pip_v7_2.dat ../../corrfunc/multipole_pair_counts/data/eBOSS_LRG_NGC_pip_v7_2.rand ../../corrfunc/multipole_pair_counts/data/eBOSS_LRG_SGC_pip_v7_2.dat ../../corrfunc/multipole_pair_counts/data/eBOSS_LRG_SGC_pip_v7_2.rand ../output/plot_data/ eBOSS_LRG_pip_v7_2 200 True

deactivate
