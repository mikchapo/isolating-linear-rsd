for i in 00-0
  do
  cp AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${i}_products/AbacusCosmos_1100box_planck_${i}_FoF_halos/z0.300/* ../abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${i}_products/AbacusCosmos_1100box_planck_${i}_FoF_halos/z0.300/
  cd ../abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${i}_products/AbacusCosmos_1100box_planck_${i}_FoF_halos/z0.300/
  tar -xzf field_subsample_ids.tar.gz
  tar -xzf field_subsamples.tar.gz
  tar -xzf halo_subsample_ids.tar.gz
  tar -xzf halo_subsamples.tar.gz
  rm *.tar.gz

  cd ~/scratch/abacus_new

  cp AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${i}_products/AbacusCosmos_1100box_planck_${i}_FoF_halos/z1.500/* ../abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${i}_products/AbacusCosmos_1100box_planck_${i}_FoF_halos/z1.500/
  cd ../abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${i}_products/AbacusCosmos_1100box_planck_${i}_FoF_halos/z1.500/
  tar -xzf field_subsample_ids.tar.gz
  tar -xzf field_subsamples.tar.gz
  tar -xzf halo_subsample_ids.tar.gz
  tar -xzf halo_subsamples.tar.gz
  rm *.tar.gz

  done
