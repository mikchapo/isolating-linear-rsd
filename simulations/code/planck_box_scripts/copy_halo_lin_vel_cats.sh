for i in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
  do
  dest_dir=~/projects/rrg-wperciva/abacus_share/AbacusCosmos_1100box_planck_rockstar_halos_z0.700_smoothed_lin_vel/
  cp ~/scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${i}_products/AbacusCosmos_1100box_planck_${i}_rockstar_halos/z0.700/halo_lin_vel_1100_tophat_5.0_smoothed_v3.dat ${dest_dir}/halo_smoothed_lin_vel_box_${i}_v3.dat
  done
