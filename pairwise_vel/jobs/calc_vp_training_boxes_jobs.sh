# Calculate the pairwise velocity distribution for the original, scaled linear, and scaled non-linear sample catalogue for the emulator
for BOXID in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
  do
  for FILL_TYPE in v-tot v-lin v-nl
    do
    export FILL_TYPE
    export INPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${BOXID}_products/AbacusCosmos_1100box_${BOXID}_rockstar_halos/z0.700/halo_lin_vel.dat"
    export OUTPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${BOXID}_products/AbacusCosmos_1100box_${BOXID}_rockstar_halos/z0.700/halo_FILL_TYPE_MASS_vp_v2.dat"
    sed "s/RUN_TIME/03:00:00/" calc_vp_template.sh > calc_vp.sh
    sed -i "s/NCPUS/32/" calc_vp.sh
    sed -i "s/JOB_NAME/calc-vp-AbacusCosmos_1100box-${BOXID}-mass/" calc_vp.sh
    sed -i "s/OUTPUT_DIR/abacus/" calc_vp.sh
    sed -i "/mail/d" calc_vp.sh
    sbatch calc_vp.sh
    done
  done
