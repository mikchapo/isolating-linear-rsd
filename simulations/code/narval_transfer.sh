# for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 planck
for i in 21
  do
  cd AbacusCosmos_1100box_${i}_products
  echo
  echo "Starting box ${i}"
  echo
  cd ic_ngc_nplt_z49.0
  scp * mj3chapm@narval.computecanada.ca:scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${i}_products/ic_ngc_nplt_z49.0/
  echo "Finished ICs"
  echo
  cd ../AbacusCosmos_1100box_${i}_FoF_halos
  scp -r * mj3chapm@narval.computecanada.ca:scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${i}_products/AbacusCosmos_1100box_${i}_FoF_halos/
  echo "Finished FoF"
  echo
  cd ../AbacusCosmos_1100box_${i}_rockstar_halos
  scp -r * mj3chapm@narval.computecanada.ca:scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${i}_products/AbacusCosmos_1100box_${i}_rockstar_halos/
  echo "Finished Rockstar"
  echo
  cd ../../
  done
