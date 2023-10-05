export INPUT_ROOT="../output/data_products/jk_pc/"
export OUTPUT_ROOT=$INPUT_ROOT
export NREG=200
export DATA_INPUT_NGC="../data/eBOSS_LRG_pip_v7_2_NGC_jk_${NREG}.dat"
export BW_INPUT_NGC="../data/eBOSS_LRG_NGC_pip_v7_2_bw.dat"
export DATA_INPUT_SGC="../data/eBOSS_LRG_pip_v7_2_SGC_jk_${NREG}.dat"
export BW_INPUT_SGC="../data/eBOSS_LRG_SGC_pip_v7_2_bw.dat"
export DD_CROSS_NORM_OUTPUT="/DD_cross_norm_eBOSS_LRG_v7_2_pip.dat"
sbatch pip/norm_count_DD_cross.sh
sleep 1
for cap in "NGC" "SGC"
  do
  export DATA_INPUT="../data/eBOSS_LRG_pip_v7_2_${cap}_jk_${NREG}.dat"
  export BW_INPUT="../data/eBOSS_LRG_${cap}_pip_v7_2_bw.dat"
  export DD_AUW_INPUT="/DD_ang_weights_eBOSS_LRG_${cap}_v7_2.dat"
  export DR_AUW_INPUT="/DR_ang_weights_eBOSS_LRG_${cap}_v7_2.dat"
  export RAND_INPUT="../data/eBOSS_LRG_pip_v7_2_${cap}_jk_${NREG}.rand"
  export DD_NORM_OUTPUT="/DD_norm_eBOSS_LRG_${cap}_v7_2_pip.dat"
  export DR_NORM_OUTPUT="/DR_norm_eBOSS_LRG_${cap}_v7_2_pip.dat"
  export RR_NORM_OUTPUT="/RR_norm_eBOSS_LRG_${cap}_v7_2_pip.dat"
#  sbatch pip/norm_count_DD.sh
#  sleep 1
#  sbatch pip/norm_count_DR.sh
#  sleep 1
#  sbatch pip/norm_count_RR.sh
#  sleep 1
  for cf in "mps" "wp"
    do
    export DD_OUTPUT="/DD_eBOSS_LRG_${cap}_v7_2_${cf}_pip_ang.dat"
    export DR_OUTPUT="/DR_eBOSS_LRG_${cap}_v7_2_${cf}_pip_ang.dat"
    export RR_OUTPUT="/RR_eBOSS_LRG_${cap}_v7_2_${cf}_pip.dat"
#    sbatch pip/DD_${cf}.sh
#    sleep 1
#    sbatch pip/DR_${cf}.sh
#    sleep 1
#    sbatch pip/RR_${cf}.sh
    done
  done


