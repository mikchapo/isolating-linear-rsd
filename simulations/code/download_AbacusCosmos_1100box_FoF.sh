wget -r -nH -np --cut-dirs=2 -R 'index.html*' --accept-regex='(_products/$|info/|_FoF_halos/$|z0.700/)'  https://nbody.rc.fas.harvard.edu/public/AbacusCosmos/AbacusCosmos_1100box_planck_products/ --no-check-certificate
find AbacusCosmos_1100box_planck_products -type f -name '*.tar.gz' -execdir tar -xzf {} \;
