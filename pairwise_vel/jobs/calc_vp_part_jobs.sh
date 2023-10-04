for BOXID in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
  do
  INPUT=/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${BOXID}_products/AbacusCosmos_1100box_${BOXID}_rockstar_halos/z0.700/ic_part_100000.dat
  OUTPUT=/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${BOXID}_products/AbacusCosmos_1100box_${BOXID}_rockstar_halos/z0.700/ic_part_100000_vp.dat
  sed "s/BOXID/${BOXID}/" calc_vp_part_template.sh > calc_vp_part.sh
  sbatch calc_vp.sh
  done
