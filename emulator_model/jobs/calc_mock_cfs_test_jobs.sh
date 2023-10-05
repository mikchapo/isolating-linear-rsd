simname=AbacusCosmos_1100box_planck
halotype=rockstar
scaling_method=split
load_dvs=False
smoothed_vel=True
subhalos=True
export simname
export boxid
export halotype
export scaling_method
export smoothed_vel
export subhalos

for boxid in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
  do
  export boxid

  smooth_type=tophat
  R_smooth=5
  export smooth_type
  export R_smooth

  for N_grid in 550 1375
    do
    export N_grid
    for scale_variable in gamma_l gamma_n
      do
      export scale_variable
      sed "s/SIMNAME/${simname}/" calc_mock_cfs_test_template.sh > calc_mock_cfs_test.sh
      sed -i "s/BOXID/${boxid}/" calc_mock_cfs_test.sh
      sed -i "s/HALOTYPE/${halotype}/" calc_mock_cfs_test.sh
      sed -i "s/SCALING_METHOD/${scaling_method}/" calc_mock_cfs_test.sh
      sed -i "s/SCALE_VARIABLE/${scale_variable}/" calc_mock_cfs_test.sh
      sed -i "s/SMOOTHED_VEL/${smoothed_vel}/" calc_mock_cfs_test.sh
      sed -i "s/N_GRID/${N_grid}/" calc_mock_cfs_test.sh
      sed -i "s/SMOOTH_TYPE/${smooth_type}/" calc_mock_cfs_test.sh
      sed -i "s/R_SMOOTH/${R_smooth}/" calc_mock_cfs_test.sh
#      sbatch calc_mock_cfs_test.sh
      done
    done

  N_grid=1100
  export N_grid

  smoothed_vel=False
  export smoothed_vel

  for scale_variable in gamma_l gamma_n
    do
    export scale_variable
    sed "s/SIMNAME/${simname}/" calc_mock_cfs_test_template.sh > calc_mock_cfs_test.sh
    sed -i "s/BOXID/${boxid}/" calc_mock_cfs_test.sh
    sed -i "s/HALOTYPE/${halotype}/" calc_mock_cfs_test.sh
    sed -i "s/SCALING_METHOD/${scaling_method}/" calc_mock_cfs_test.sh
    sed -i "s/SCALE_VARIABLE/${scale_variable}/" calc_mock_cfs_test.sh
    sed -i "s/SMOOTHED_VEL/${smoothed_vel}/" calc_mock_cfs_test.sh
    sed -i "s/N_GRID/${N_grid}/" calc_mock_cfs_test.sh
    sed -i "s/SMOOTH_TYPE/${smooth_type}/" calc_mock_cfs_test.sh
    sed -i "s/R_SMOOTH/${R_smooth}/" calc_mock_cfs_test.sh
    sbatch calc_mock_cfs_test.sh
    done

  smoothed_vel=True
  export smoothed_vel

  for R_smooth in 3 5 7
    do
    export R_smooth
    for scale_variable in gamma_l gamma_n
      do
      export scale_variable
      sed "s/SIMNAME/${simname}/" calc_mock_cfs_test_template.sh > calc_mock_cfs_test.sh
      sed -i "s/BOXID/${boxid}/" calc_mock_cfs_test.sh
      sed -i "s/HALOTYPE/${halotype}/" calc_mock_cfs_test.sh
      sed -i "s/SCALING_METHOD/${scaling_method}/" calc_mock_cfs_test.sh
      sed -i "s/SCALE_VARIABLE/${scale_variable}/" calc_mock_cfs_test.sh
      sed -i "s/SMOOTHED_VEL/${smoothed_vel}/" calc_mock_cfs_test.sh
      sed -i "s/N_GRID/${N_grid}/" calc_mock_cfs_test.sh
      sed -i "s/SMOOTH_TYPE/${smooth_type}/" calc_mock_cfs_test.sh
      sed -i "s/R_SMOOTH/${R_smooth}/" calc_mock_cfs_test.sh
#      sbatch calc_mock_cfs_test.sh
      done
    done

  smooth_type=gaussian
  R_smooth=2
  export smooth_type
  export R_smooth

  for scale_variable in gamma_l gamma_n
    do
    export scale_variable
    sed "s/SIMNAME/${simname}/" calc_mock_cfs_test_template.sh > calc_mock_cfs_test.sh
    sed -i "s/BOXID/${boxid}/" calc_mock_cfs_test.sh
    sed -i "s/HALOTYPE/${halotype}/" calc_mock_cfs_test.sh
    sed -i "s/SCALING_METHOD/${scaling_method}/" calc_mock_cfs_test.sh
    sed -i "s/SCALE_VARIABLE/${scale_variable}/" calc_mock_cfs_test.sh
    sed -i "s/SMOOTHED_VEL/${smoothed_vel}/" calc_mock_cfs_test.sh
    sed -i "s/N_GRID/${N_grid}/" calc_mock_cfs_test.sh
    sed -i "s/SMOOTH_TYPE/${smooth_type}/" calc_mock_cfs_test.sh
    sed -i "s/R_SMOOTH/${R_smooth}/" calc_mock_cfs_test.sh
#    sbatch calc_mock_cfs_test.sh
    done
  done
