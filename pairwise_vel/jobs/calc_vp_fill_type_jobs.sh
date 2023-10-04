# Calculate the pairwise velocity distribution from the total, linear, and non-linear velocity components of a sample Abacus box
# for IC_TYPE in default low_redshift
for IC_TYPE in ic_no_growth_corr
# for IC_TYPE in low_redshift
  do
#   for FILL_TYPE in v-tot v-lin v-nl
  for FILL_TYPE in v-tot v-lin v-nl
    do
    export INPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/halos_lin_vel_${IC_TYPE}.dat"
    export OUTPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/FoF_${IC_TYPE}_FILL_TYPE_MASS_vp.dat"
    export FILL_TYPE
    sed "s/RUN_TIME/01:30:00/" calc_vp_template.sh > calc_vp.sh
    sed -i "s/NCPUS/32/" calc_vp.sh
    sed -i "s/JOB_NAME/calc-vp-abacus-planck-z0.70-${IC_TYPE}-${FILL_TYPE}-mass/" calc_vp.sh
    sed -i "s/OUTPUT_DIR/abacus/" calc_vp.sh
    sbatch calc_vp.sh
#    echo $INPUT $OUTPUT $FILL_TYPE| ../code/01-calc_vp_halo_llist_pbc.exe
    done
  done
