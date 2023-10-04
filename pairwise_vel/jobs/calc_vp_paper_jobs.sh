# Calculate the pairwise velocity distribution for the original, scaled linear, and scaled non-linear sample catalogue for the emulator
for BOXID in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
  do
  for FILL_TYPE in v-lin v-nl
    do
    export FILL_TYPE
    export INPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${BOXID}_products/AbacusCosmos_1100box_planck_${BOXID}_rockstar_halos/z0.700/halo_lin_vel_all-halos.dat"
    export OUTPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${BOXID}_products/AbacusCosmos_1100box_planck_${BOXID}_rockstar_halos/z0.700/pairwise_vel/halo_FILL_TYPE_MASS_vp_v2_all-halos.dat"
    sed "s/RUN_TIME/03:00:00/" calc_vp_template.sh > calc_vp.sh
    sed -i "s/NCPUS/32/" calc_vp.sh
    sed -i "s/JOB_NAME/calc-vp-test-${BOXID}-${FILL_TYPE}-mass/" calc_vp.sh
    sed -i "s/OUTPUT_DIR/abacus/" calc_vp.sh
    sbatch calc_vp.sh
    done
  done

for BOXID in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
  do
  for FILL_TYPE in v-lin v-nl
    do
    export FILL_TYPE
    export INPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${BOXID}_products/AbacusCosmos_1100box_planck_${BOXID}_rockstar_halos/z0.700/halo_lin_vel_1100_tophat_5.0_smoothed_v2_all-halos.dat"
    export OUTPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${BOXID}_products/AbacusCosmos_1100box_planck_${BOXID}_rockstar_halos/z0.700/pairwise_vel/halo_FILL_TYPE_MASS_vp_v2_vel-sm_all-halos.dat"
    sed "s/RUN_TIME/03:00:00/" calc_vp_template.sh > calc_vp.sh
    sed -i "s/NCPUS/32/" calc_vp.sh
    sed -i "s/JOB_NAME/calc-vp-test-${BOXID}-${FILL_TYPE}-mass-vel-sm/" calc_vp.sh
    sed -i "s/OUTPUT_DIR/abacus/" calc_vp.sh
    sbatch calc_vp.sh
    done
  done

for BOXID in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
  do
  FILL_TYPE=v-tot
  export FILL_TYPE
  export INPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${BOXID}_products/AbacusCosmos_1100box_planck_${BOXID}_rockstar_halos/z0.700/halo_lin_vel_all-halos.dat"
  export OUTPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${BOXID}_products/AbacusCosmos_1100box_planck_${BOXID}_rockstar_halos/z0.700/pairwise_vel/halo_FILL_TYPE_MASS_vp_v2_all-halos.dat"
  sed "s/RUN_TIME/03:00:00/" calc_vp_template.sh > calc_vp.sh
  sed -i "s/NCPUS/32/" calc_vp.sh
  sed -i "s/JOB_NAME/calc-vp-test-${BOXID}-${FILL_TYPE}-mass/" calc_vp.sh
  sed -i "s/OUTPUT_DIR/abacus/" calc_vp.sh
  sbatch calc_vp.sh
  done
