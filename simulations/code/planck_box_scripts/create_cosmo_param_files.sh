for i in {0..19}
  do
  cd ~/scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_00-${i}_products/AbacusCosmos_1100box_planck_00-${i}_rockstar_halos/info/
  sed -n '/^H0 = / s/.*\= *//p' abacus.par > cosmo_params.dat
  sed -n '/^omch2 = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/^ombh2 = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/^ns = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/^sigma_8 = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/^N_eff = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/^w0 = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/^ZD_Pk_sigma = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  done
