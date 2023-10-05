export OUTPUT_ROOT="../output/data_products/jk_pc/"
export NREG=200
for cap in "NGC" "SGC"
  do
  export DATA_INPUT="../data/eBOSS_LRG_pip_v7_2_${cap}_jk_200.dat"
  export BW_INPUT="../data/eBOSS_LRG_${cap}_pip_v7_2_bw.dat"
  export RAND_INPUT="../data/eBOSS_LRG_pip_v7_2_${cap}_jk_200.rand"
#  for pc in "DD" "DR"
  for pc in "DD"
    do
    export NORM_OUTPUT="/${pc}_norm_eBOSS_LRG_${cap}_v7_2_ang.dat"
    sbatch ang/norm_count_${pc}.sh
    sleep 1
    export PAR_OUTPUT="/${pc}_eBOSS_LRG_${cap}_v7_2_ang_par.dat"
    export FIB_OUTPUT="/${pc}_eBOSS_LRG_${cap}_v7_2_ang_fib_pip.dat"
#    sbatch ang/${pc}.sh
    sleep 1
    done
  done


