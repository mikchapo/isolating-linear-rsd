simname=AbacusCosmos_1100box
boxid=00
halotype=rockstar
scaling_method=split
load_dvs=False
smoothed_vel=True
export simname
export boxid
export halotype
export scaling_method
export load_dvs
export smoothed_vel

smooth_type=tophat
R_smooth=5
export smooth_type
export R_smooth

for N_grid in 550 1100 1375
  do
  export N_grid
  for scale_variable in gamma_l gamma_n
    do
    export scale_variable
    sed "s/SIMNAME/${simname}/" mock_cfs_template.sh > mock_cfs.sh
    sed -i "s/BOXID/${boxid}/" mock_cfs.sh
    sed -i "s/HALOTYPE/${halotype}/" mock_cfs.sh
    sed -i "s/SCALING_METHOD/${scaling_method}/" mock_cfs.sh
    sed -i "s/SCALE_VARIABLE/${scale_variable}/" mock_cfs.sh
    sed -i "s/SMOOTHED_VEL/${smoothed_vel}/" mock_cfs.sh
    sed -i "s/N_GRID/${N_grid}/" mock_cfs.sh
    sed -i "s/SMOOTH_TYPE/${smooth_type}/" mock_cfs.sh
    sed -i "s/R_SMOOTH/${R_smooth}/" mock_cfs.sh
#    sbatch mock_cfs.sh
    done
  done

N_grid=1100
export N_grid

#for R_smooth in 1 2 3 4 6
for R_smooth in 1 2 3 4
  do
  export R_smooth
  for scale_variable in gamma_l gamma_n
    do
    export scale_variable
    sed "s/SIMNAME/${simname}/" mock_cfs_template.sh > mock_cfs.sh
    sed -i "s/BOXID/${boxid}/" mock_cfs.sh
    sed -i "s/HALOTYPE/${halotype}/" mock_cfs.sh
    sed -i "s/SCALING_METHOD/${scaling_method}/" mock_cfs.sh
    sed -i "s/SCALE_VARIABLE/${scale_variable}/" mock_cfs.sh
    sed -i "s/SMOOTHED_VEL/${smoothed_vel}/" mock_cfs.sh
    sed -i "s/N_GRID/${N_grid}/" mock_cfs.sh
    sed -i "s/SMOOTH_TYPE/${smooth_type}/" mock_cfs.sh
    sed -i "s/R_SMOOTH/${R_smooth}/" mock_cfs.sh
    sbatch mock_cfs.sh
    done
  done

smooth_type=gaussian
export smooth_type

for R_smooth in 1 2
  do
  export R_smooth
  for scale_variable in gamma_l gamma_n
    do
    export scale_variable
    sed "s/SIMNAME/${simname}/" mock_cfs_template.sh > mock_cfs.sh
    sed -i "s/BOXID/${boxid}/" mock_cfs.sh
    sed -i "s/HALOTYPE/${halotype}/" mock_cfs.sh
    sed -i "s/SCALING_METHOD/${scaling_method}/" mock_cfs.sh
    sed -i "s/SCALE_VARIABLE/${scale_variable}/" mock_cfs.sh
    sed -i "s/SMOOTHED_VEL/${smoothed_vel}/" mock_cfs.sh
    sed -i "s/N_GRID/${N_grid}/" mock_cfs.sh
    sed -i "s/SMOOTH_TYPE/${smooth_type}/" mock_cfs.sh
    sed -i "s/R_SMOOTH/${R_smooth}/" mock_cfs.sh
#    sbatch mock_cfs.sh
    done
  done
