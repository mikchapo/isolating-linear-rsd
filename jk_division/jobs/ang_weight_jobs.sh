export INPUT_ROOT="../output/data_products/jk_pc/"
export OUTPUT_ROOT=$INPUT_ROOT
export NREG=200
for cap in "NGC" "SGC"
  do
  for pc in "DD" "DR"
    do
    export PAR_INPUT="/${pc}_eBOSS_LRG_${cap}_v7_2_ang_par.dat"
    export FIB_INPUT="/${pc}_eBOSS_LRG_${cap}_v7_2_ang_fib_pip.dat"
    export NORM_INPUT="/${pc}_norm_eBOSS_LRG_${cap}_v7_2_ang.dat"
    export WEIGHT_OUTPUT="/${pc}_ang_weights_eBOSS_LRG_${cap}_v7_2.dat"
    sbatch calc_ang_weight.sh
    sleep 1
    done
  done


