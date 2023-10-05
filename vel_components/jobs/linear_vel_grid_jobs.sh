simname=AbacusCosmos_1100box
N_grid=1100
export simname
export N_grid

for boxid in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 planck
  do
  export boxid
  sed "s/SIMNAME/${simname}/" linear_vel_grid_template.sh > linear_vel_grid.sh
  sed -i "s/MEMORY/70G/" linear_vel_grid.sh
  sed -i "s/BOXID/${boxid}/" linear_vel_grid.sh
  sed -i "s/N_GRID/${N_grid}/" linear_vel_grid.sh
  sed -i "s/OUTPUT_DIR/training/" linear_vel_grid.sh
  sbatch linear_vel_grid.sh
  done

simname=AbacusCosmos_1100box_planck
export simname

for boxid in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
  do
  export boxid
  sed "s/SIMNAME/${simname}/" linear_vel_grid_template.sh > linear_vel_grid.sh
  sed -i "s/MEMORY/70G/" linear_vel_grid.sh
  sed -i "s/BOXID/${boxid}/" linear_vel_grid.sh
  sed -i "s/N_GRID/${N_grid}/" linear_vel_grid.sh
  sed -i "s/OUTPUT_DIR/test/" linear_vel_grid.sh
  sbatch linear_vel_grid.sh
  done
