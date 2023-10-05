# for i in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
for i in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
  do
#   mkdir AbacusCosmos_1100box_planck_${i}_products
#   cd AbacusCosmos_1100box_planck_${i}_products
#   mkdir ic_ngc_nplt_z49.0
#   mkdir AbacusCosmos_1100box_planck_${i}_FoF_halos
#   mkdir AbacusCosmos_1100box_planck_${i}_rockstar_halos

  cd AbacusCosmos_1100box_planck_${i}_products/AbacusCosmos_1100box_planck_${i}_rockstar_halos/z0.700/
  echo $i
  ls pairwise_vel/
  cd ../../../

#  cd AbacusCosmos_1100box_planck_${i}_products/AbacusCosmos_1100box_planck_${i}_rockstar_halos/z0.700/pairwise_vel
#  mkdir v1
#  rm *.dat
#  cd ../../../../

  done
