# Calculate the pairwise velocity distribution from the total, linear, and non-linear velocity components of a sample Abacus box
for IC_TYPE in default no_growth_corr ngc_nplt ngc_nplt_0.7
# for IC_TYPE in low_redshift
  do
  for FILL_TYPE in v-tot v-lin v-nl
    do
    export INPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/halo_lin_vel_${IC_TYPE}.dat"
    export OUTPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/FoF_${IC_TYPE}_test_FILL_TYPE_MASS_vp.dat"
    export FILL_TYPE
    sed "s/RUN_TIME/01:30:00/" calc_vp_template.sh > calc_vp.sh
    sed -i "s/NCPUS/32/" calc_vp.sh
    sed -i "s/JOB_NAME/calc-vp-abacus-planck-z0.70-${IC_TYPE}-${FILL_TYPE}-mass/" calc_vp.sh
    sed -i "s/OUTPUT_DIR/abacus/" calc_vp.sh
    sbatch calc_vp.sh
    done
  done
