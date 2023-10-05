for i in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
# for i in 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
  do
  cd AbacusCosmos_1100box_planck_${i}_products
  echo
  echo "Starting box ${i}"
#  echo
#  cd ic_ngc_nplt_z49.0
#  scp * mj3chapm@narval.computecanada.ca:scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${i}_products/ic_ngc_nplt_z49.0/
#  echo "Finished ICs"
  echo
  cd AbacusCosmos_1100box_planck_${i}_FoF_halos
  scp -r * mj3chapm@narval.computecanada.ca:scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${i}_products/AbacusCosmos_1100box_planck_${i}_FoF_halos/
  echo "Finished FoF"
#   echo
#   cd AbacusCosmos_1100box_planck_${i}_rockstar_halos
#   scp -r * mj3chapm@narval.computecanada.ca:scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${i}_products/AbacusCosmos_1100box_planck_${i}_rockstar_halos/
#   echo "Finished Rockstar"
  echo
  cd ../../
  done
