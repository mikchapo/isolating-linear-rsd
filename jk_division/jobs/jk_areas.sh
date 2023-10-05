#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=jk-area-pip
#SBATCH --output=../output/%x-%j.out

module load python/3.6

source ~/env/bin/activate

python ../code/jk_areas.py ../data/eBOSS_LRG_NGC_pip_v7_2.dat ../data/eBOSS_LRG_NGC_pip_v7_2.rand ../data/eBOSS_LRG_SGC_pip_v7_2.dat ../data/eBOSS_LRG_SGC_pip_v7_2.rand ../output/eBOSS_LRG_pip_v7_2 200

deactivate
