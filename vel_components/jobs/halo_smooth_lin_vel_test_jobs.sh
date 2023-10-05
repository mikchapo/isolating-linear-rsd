simname=AbacusCosmos_1100box_planck
halo_type=rockstar
subhalos=False
export simname
export halo_type
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
    sed "s/SIMNAME/${simname}/" halo_smooth_lin_vel_template.sh > halo_smooth_lin_vel.sh
    sed -i "s/BOXID/${boxid}/" halo_smooth_lin_vel.sh
    sed -i "s/HALO_TYPE/${halo_type}/" halo_smooth_lin_vel.sh
    sed -i "s/N_GRID/${N_grid}/" halo_smooth_lin_vel.sh
    sed -i "s/SMOOTH_TYPE/${smooth_type}/" halo_smooth_lin_vel.sh
    sed -i "s/R_SMOOTH/${R_smooth}/" halo_smooth_lin_vel.sh
    sed -i "s/OUTPUT_DIR/test/" halo_smooth_lin_vel.sh
    sbatch halo_smooth_lin_vel.sh
    done

  N_grid=1100
  export N_grid

  for R_smooth in 3 7
    do
    export R_smooth
    sed "s/SIMNAME/${simname}/" halo_smooth_lin_vel_template.sh > halo_smooth_lin_vel.sh
    sed -i "s/BOXID/${boxid}/" halo_smooth_lin_vel.sh
    sed -i "s/HALO_TYPE/${halo_type}/" halo_smooth_lin_vel.sh
    sed -i "s/N_GRID/${N_grid}/" halo_smooth_lin_vel.sh
    sed -i "s/SMOOTH_TYPE/${smooth_type}/" halo_smooth_lin_vel.sh
    sed -i "s/R_SMOOTH/${R_smooth}/" halo_smooth_lin_vel.sh
    sed -i "s/OUTPUT_DIR/test/" halo_smooth_lin_vel.sh
    sbatch halo_smooth_lin_vel.sh
    done

  smooth_type=gaussian
  R_smooth=2
  export smooth_type
  export R_smooth

  sed "s/SIMNAME/${simname}/" halo_smooth_lin_vel_template.sh > halo_smooth_lin_vel.sh
  sed -i "s/BOXID/${boxid}/" halo_smooth_lin_vel.sh
  sed -i "s/HALO_TYPE/${halo_type}/" halo_smooth_lin_vel.sh
  sed -i "s/N_GRID/${N_grid}/" halo_smooth_lin_vel.sh
  sed -i "s/SMOOTH_TYPE/${smooth_type}/" halo_smooth_lin_vel.sh
  sed -i "s/R_SMOOTH/${R_smooth}/" halo_smooth_lin_vel.sh
  sed -i "s/OUTPUT_DIR/test/" halo_smooth_lin_vel.sh
  sbatch halo_smooth_lin_vel.sh

  smooth_type=tophat
  R_smooth=5
  subhalos=True
  export smooth_type
  export R_smooth
  export subhalos

  sed "s/SIMNAME/${simname}/" halo_smooth_lin_vel_template.sh > halo_smooth_lin_vel.sh
  sed -i "s/BOXID/${boxid}/" halo_smooth_lin_vel.sh
  sed -i "s/HALO_TYPE/${halo_type}/" halo_smooth_lin_vel.sh
  sed -i "s/N_GRID/${N_grid}/" halo_smooth_lin_vel.sh
  sed -i "s/SMOOTH_TYPE/${smooth_type}/" halo_smooth_lin_vel.sh
  sed -i "s/R_SMOOTH/${R_smooth}/" halo_smooth_lin_vel.sh
  sed -i "s/OUTPUT_DIR/test/" halo_smooth_lin_vel.sh
  sbatch halo_smooth_lin_vel.sh

  done
