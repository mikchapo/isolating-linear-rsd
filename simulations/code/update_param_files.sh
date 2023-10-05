for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 planck
  do
  cd ~/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${i}_products/AbacusCosmos_1100box_${i}_rockstar_halos/info/
  # Use % instead of / to divide sed command so that a path containing / can be replaced
  sed "s%/localscratch/lhgarrison/AbacusCosmos_1100box_${i}/ic%/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${i}_products/ic_ngc_nplt_z49.0%" abacus.par > abacus_ngc_nplt.par
  sed -i "s%/dev/shm/lhgarrison/abacus//zeldovich-PLT/eigmodes128%/home/mj3chapm/P2/vel_components/zeldovich-PLT/eigmodes128%" abacus_ngc_nplt.par
  sed -i "s%info/camb_matterpower.dat%/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${i}_products/AbacusCosmos_1100box_${i}_rockstar_halos/info/camb_matterpower.dat%" abacus_ngc_nplt.par
  sed -i '/ZD_qPLT_rescale = 1/a\ZD_Version = 1' abacus_ngc_nplt.par
  sed -i "s%ZD_qPLT = 1%ZD_qPLT = 0%" abacus_ngc_nplt.par
  sed -i "s%ZD_qPLT_rescale = 1%ZD_qPLT_rescale = 0%" abacus_ngc_nplt.par
  done

i=planck
cd ~/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${i}_products/AbacusCosmos_1100box_${i}_rockstar_halos/info/
# Use % instead of / to divide sed command so that a path containing / can be replaced
sed "s%/mnt/raid/lgarrison//AbacusCosmos_1100box_planck/ic%/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${i}_products/ic_ngc_nplt_z49.0%" abacus.par > abacus_ngc_nplt.par
sed -i "s%/home/lgarrison/abacus/zeldovich-PLT/eigmodes128%/home/mj3chapm/P2/vel_components/zeldovich-PLT/eigmodes128%" abacus_ngc_nplt.par
sed -i "s%info/camb_matterpower.dat%/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${i}_products/AbacusCosmos_1100box_${i}_rockstar_halos/info/camb_matterpower.dat%" abacus_ngc_nplt.par
sed -i '/ZD_qPLT_rescale = 1/a\ZD_Version = 1' abacus_ngc_nplt.par
sed -i "s%ZD_qPLT = 1%ZD_qPLT = 0%" abacus_ngc_nplt.par
sed -i "s%ZD_qPLT_rescale = 1%ZD_qPLT_rescale = 0%" abacus_ngc_nplt.par
