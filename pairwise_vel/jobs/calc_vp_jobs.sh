for BOX_ID in 000 001 002 003 004 005 006
# for BOX_ID in 000
  do
for REALIZATION in 000 001 002 003 004
#  for REALIZATION in 000
    do
    export INPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/halos/aemulus_z0.70_TestBox${BOX_ID}-${REALIZATION}_m200b_halos.list"
    export OUTPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_TestBox${BOX_ID}-${REALIZATION}/m200b_MASS_vp.dat"
    sed "s/RUN_TIME/01:00:00/" calc_vp_template.sh > calc_vp.sh
    sed -i "s/NCPUS/32/" calc_vp.sh
    sed -i "s/JOB_NAME/calc-vp-TestBox${BOX_ID}-${REALIZATION}-z0.70-mass/" calc_vp.sh
    sbatch calc_vp.sh
    done
  done
