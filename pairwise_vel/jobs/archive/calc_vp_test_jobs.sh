# Testing tyhe pairwise velocity code to see where it is going wrong
export FILL_TYPE=v-nl
export INPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/halos_lin_vel.dat"
# for SCRIPT_NAME in 01-calc_vp_halo_llist_pbc.exe 02-calc_vp_halo_llist_pbc_final-var.exe 03-calc_vp_halo_llist_pbc_dyn-var.exe
for SCRIPT_NAME in 01-calc_vp_halo_llist_pbc.exe 02-calc_vp_halo_llist_pbc_final-var.exe
  do
#  for NBIN in 40 80
#    do
  export OUTPUT="/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/FoF_FILL_TYPE_MASS_vp_script${SCRIPT_NAME:1:1}_maxs-10.dat"
#    export NBIN
  sed "s/RUN_TIME/00:30:00/" calc_vp_test_template.sh > calc_vp_test.sh
  sed -i "s/NCPUS/32/" calc_vp_test.sh
  sed -i "s/JOB_NAME/calc-vp-abacus-planck-z0.70-${FILL_TYPE}-mass-script${SCRIPT_NAME:1:1}-maxs-10/" calc_vp_test.sh
  sed -i "s/OUTPUT_DIR/abacus/" calc_vp_test.sh
  sed -i "s/SCRIPT/${SCRIPT_NAME}/" calc_vp_test.sh
  sbatch calc_vp_test.sh
#    done
  done
