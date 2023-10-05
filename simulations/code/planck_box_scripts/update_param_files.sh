for i in {0..19}
  do
  cd ~/scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_00-${i}_products/AbacusCosmos_1100box_planck_00-${i}_rockstar_halos/info/
  # Use % instead of / to divide sed command so that a path containing / can be replaced
  sed "s%/localscratch/lhgarrison/AbacusCosmos_1100box_planck_00-${i}/ic%/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_00-${i}_products/ic_ngc_nplt_z49.0%" abacus.par > abacus_ngc_nplt.par
  sed -i "s%/dev/shm/lhgarrison/abacus//zeldovich-PLT/eigmodes128%/home/mj3chapm/P2/vel_components/zeldovich-PLT/eigmodes128%" abacus_ngc_nplt.par
  sed -i "s%info/camb_matterpower.dat%/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_00-${i}_products/AbacusCosmos_1100box_planck_00-${i}_rockstar_halos/info/camb_matterpower.dat%" abacus_ngc_nplt.par
  sed -i '/ZD_qPLT_rescale = 1/a\ZD_Version = 1' abacus_ngc_nplt.par
  sed -i "s%ZD_qPLT = 1%ZD_qPLT = 0%" abacus_ngc_nplt.par
  sed -i "s%ZD_qPLT_rescale = 1%ZD_qPLT_rescale = 0%" abacus_ngc_nplt.par
  done
