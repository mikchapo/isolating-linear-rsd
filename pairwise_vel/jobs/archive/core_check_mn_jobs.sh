REALIZATION=000
export INPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/halos/aemulus_z0.70_TestBox001-${REALIZATION}_m200b_halos.list"
export M=100
export OUTPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_TestBox001-${REALIZATION}/m200b_MASS_vp.dat"
sed "s/RUN_TIME/01:00:00/" calc_vp_mn_template.sh > calc_vp_mn.sh
sed -i "s/NCPUS/32/" calc_vp_mn.sh
sed -i "s/JOB_NAME/calc-vp-TestBox001-${REALIZATION}-z0.70-mass/" calc_vp_mn.sh
sbatch calc_vp_mn.sh


REALIZATION=001
export INPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/halos/aemulus_z0.70_TestBox001-${REALIZATION}_m200b_halos.list"
export M=100
export OUTPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_TestBox001-${REALIZATION}/m200b_MASS_vp.dat"
sed "s/RUN_TIME/02:00:00/" calc_vp_mn_template.sh > calc_vp_mn.sh
sed -i "s/NCPUS/16/" calc_vp_mn.sh
sed -i "s/JOB_NAME/calc-vp-TestBox001-${REALIZATION}-z0.70-mass/" calc_vp_mn.sh
sbatch calc_vp_mn.sh


REALIZATION=002
export INPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/halos/aemulus_z0.70_TestBox001-${REALIZATION}_m200b_halos.list"
export M=100
export OUTPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_TestBox001-${REALIZATION}/m200b_MASS_vp.dat"
sed "s/RUN_TIME/04:00:00/" calc_vp_mn_template.sh > calc_vp_mn.sh
sed -i "s/NCPUS/8/" calc_vp_mn.sh
sed -i "s/JOB_NAME/calc-vp-TestBox001-${REALIZATION}-z0.70-mass/" calc_vp_mn.sh
sbatch calc_vp_mn.sh


REALIZATION=003
export INPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/halos/aemulus_z0.70_TestBox001-${REALIZATION}_m200b_halos.list"
export M=100
export OUTPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_TestBox001-${REALIZATION}/m200b_MASS_vp.dat"
sed "s/RUN_TIME/08:00:00/" calc_vp_mn_template.sh > calc_vp_mn.sh
sed -i "s/NCPUS/4/" calc_vp_mn.sh
sed -i "s/JOB_NAME/calc-vp-TestBox001-${REALIZATION}-z0.70-mass/" calc_vp_mn.sh
sbatch calc_vp_mn.sh


REALIZATION=004
export INPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/halos/aemulus_z0.70_TestBox001-${REALIZATION}_m200b_halos.list"
export M=100
export OUTPUT="../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_TestBox001-${REALIZATION}/m200b_MASS_vp.dat"
sed "s/RUN_TIME/16:00:00/" calc_vp_mn_template.sh > calc_vp_mn.sh
sed -i "s/NCPUS/2/" calc_vp_mn.sh
sed -i "s/JOB_NAME/calc-vp-TestBox001-${REALIZATION}-z0.70-mass/" calc_vp_mn.sh
sbatch calc_vp_mn.sh


