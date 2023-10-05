for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 planck
  do
  cd ~/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${i}_products/AbacusCosmos_1100box_${i}_rockstar_halos/info/
  sed -n '/H0 = / s/.*\= *//p' abacus.par > cosmo_params.dat
  sed -n '/omch2 = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/ombh2 = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/ns = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/sigma_8 = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/N_eff = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/w0 = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  sed -n '/ZD_Pk_sigma = / s/.*\= *//p' abacus.par >> cosmo_params.dat
  done
