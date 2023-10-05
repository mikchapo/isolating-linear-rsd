simname=AbacusCosmos_1100box_planck
halotype=rockstar
scaling_method=split
smoothed_vel=True
subhalos=True
export simname
export boxid
export halotype
export scaling_method
export smoothed_vel
export subhalos

boxid=00-2
export boxid

N_grid=1100
smooth_type=tophat
R_smooth=5
export smooth_type
export R_smooth
export N_grid

scale_variable=gamma_n
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
#sbatch calc_mock_cfs_test.sh

R_smooth=7
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
#  sbatch calc_mock_cfs_test.sh
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
#  sbatch calc_mock_cfs_test.sh
  done



boxid=00-3
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
    sbatch calc_mock_cfs_test.sh
    done
  done

N_grid=1100
export N_grid

for R_smooth in 3 5
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
    sbatch calc_mock_cfs_test.sh
    done
  done


boxid=00-6
export boxid

smooth_type=tophat
R_smooth=7
N_grid=1100
export N_grid
export smooth_type
export R_smooth


scale_variable=gamma_n
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


smooth_type=gaussian
R_smooth=2
export smooth_type
export R_smooth

scale_variable=gamma_l
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
