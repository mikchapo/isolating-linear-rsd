#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=rrg-wperciva
#SBATCH --mail-user=mj3chapm@uwaterloo.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=download-sims-AbacusCosmos_1100box_planck_00-0-z-0.3-1.5-particles
#SBATCH --output=../output/job_logs/%x-%j.out

# cd ~/scratch/abacus
cd ~/scratch/abacus_new

# wget -r -nH -np --cut-dirs=2 -R 'index.html*' --accept-regex='(_products/$|info/|_FoF_halos/$|z0.700/)'  https://nbody.rc.fas.harvard.edu/public/AbacusCosmos/AbacusCosmos_1100box_planck_products/ --no-check-certificate
# wget -r -nH -np --cut-dirs=2 -R 'index.html*' --accept-regex='(_products/$|info/|_FoF_halos/$|z0.700/)' --reject-regex='(FoF_halos/.*/halos.tar.gz|FoF_halos/.*/halo_subsamples.tar.gz|FoF_halos/.*/halo_subsample_ids.tar.gz)' https://nbody.rc.fas.harvard.edu/public/AbacusCosmos/AbacusCosmos_1100box_products/ --no-check-certificate
wget -r -nH -np --cut-dirs=2 -R 'index.html*' --accept-regex='(_products/$|info/|_FoF_halos/$|z0.300/|z1.500/)' --reject-regex='(FoF_halos/.*/halos.tar.gz)' https://nbody.rc.fas.harvard.edu/public/AbacusCosmos/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_00-0_products/
# find AbacusCosmos_1100box_products -type f -name '*.tar.gz' -execdir tar -xzf {} \;
