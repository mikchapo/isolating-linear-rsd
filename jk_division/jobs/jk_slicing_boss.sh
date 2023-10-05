#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-wperciva
#SBATCH --job-name=jk-slicing-boss
#SBATCH --output=../output/%x-%j.out

module load python/3.6

source ~/env/bin/activate

data_file_ngc="../../../projects/rrg-wperciva/BOSS_mike/boss_x_eboss_cats/cmass_NGC_DR12_masked.fits"
rand_file_ngc="../../../projects/rrg-wperciva/BOSS_mike/boss_x_eboss_cats/cmass+eboss_NGC_rand_DR12_masked.fits"
data_file_sgc="../../../projects/rrg-wperciva/BOSS_mike/boss_x_eboss_cats/cmass_SGC_DR12_masked.fits"
rand_file_sgc="../../../projects/rrg-wperciva/BOSS_mike/boss_x_eboss_cats/cmass+eboss_SGC_rand_DR12_masked.fits"
output_root="../../../projects/rrg-wperciva/eBOSS_faizan/boss_x_eboss_cats/cmass_DR12_masked"

python ../code/jk_slicing.py  $data_file_ngc $rand_file_ngc $data_file_sgc $rand_file_sgc $output_root 200 0.59049

deactivate
