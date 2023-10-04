# Calculate the pairwise velocity distribution for the original, scaled linear, and scaled non-linear sample catalogue for the emulator
for HALO_TYPE in rockstar FoF
  do
  for FILL_TYPE in v-tot v-lin v-nl
    do
    export FILL_TYPE
    export INPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${HALO_TYPE}_halos/z0.700/halo_lin_vel.dat"
    export OUTPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_${HALO_TYPE}_halos/z0.700/halo_FILL_TYPE_MASS_vp_v2.dat"
    sed "s/RUN_TIME/03:00:00/" calc_vp_template.sh > calc_vp.sh
    sed -i "s/NCPUS/32/" calc_vp.sh
    sed -i "s/JOB_NAME/calc-vp-check-${HALO_TYPE}-mass/" calc_vp.sh
    sed -i "s/OUTPUT_DIR/abacus/" calc_vp.sh
    sbatch calc_vp.sh
    done
  done
