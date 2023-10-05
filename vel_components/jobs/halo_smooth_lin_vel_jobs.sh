simname=AbacusCosmos_1100box
N_grid=1100
halo_type=rockstar
smooth_type=tophat
R_smooth=5
subhalos=False
export simname
export N_grid
export halo_type
export smooth_type
export R_smooth
export subhalos

for boxid in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 planck
  do
  export boxid
  sed "s/SIMNAME/${simname}/" halo_smooth_lin_vel_template.sh > halo_smooth_lin_vel.sh
  sed -i "s/BOXID/${boxid}/" halo_smooth_lin_vel.sh
  sed -i "s/HALO_TYPE/${halo_type}/" halo_smooth_lin_vel.sh
  sed -i "s/N_GRID/${N_grid}/" halo_smooth_lin_vel.sh
  sed -i "s/SMOOTH_TYPE/${smooth_type}/" halo_smooth_lin_vel.sh
  sed -i "s/R_SMOOTH/${R_smooth}/" halo_smooth_lin_vel.sh
  sed -i "s/OUTPUT_DIR/training/" halo_smooth_lin_vel.sh
  sbatch halo_smooth_lin_vel.sh
  done

simname=AbacusCosmos_1100box_planck
export simname

for boxid in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
# for boxid in 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
# for boxid in 00-0
  do
  export boxid
  sed "s/SIMNAME/${simname}/" halo_smooth_lin_vel_template.sh > halo_smooth_lin_vel.sh
  sed -i "s/BOXID/${boxid}/" halo_smooth_lin_vel.sh
  sed -i "s/HALO_TYPE/${halo_type}/" halo_smooth_lin_vel.sh
  sed -i "s/N_GRID/${N_grid}/" halo_smooth_lin_vel.sh
  sed -i "s/SMOOTH_TYPE/${smooth_type}/" halo_smooth_lin_vel.sh
  sed -i "s/R_SMOOTH/${R_smooth}/" halo_smooth_lin_vel.sh
  sed -i "s/OUTPUT_DIR/test/" halo_smooth_lin_vel.sh
#  sbatch halo_smooth_lin_vel.sh
  done
