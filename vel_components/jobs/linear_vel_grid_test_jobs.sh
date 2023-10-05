simname=AbacusCosmos_1100box_planck
export simname

for boxid in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
  do
  export boxid

  N_grid=550
  export N_grid
  sed "s/SIMNAME/${simname}/" linear_vel_grid_template.sh > linear_vel_grid.sh
  sed -i "s/MEMORY/40G/" linear_vel_grid.sh
  sed -i "s/BOXID/${boxid}/" linear_vel_grid.sh
  sed -i "s/N_GRID/${N_grid}/" linear_vel_grid.sh
  sed -i "s/OUTPUT_DIR/test/" linear_vel_grid.sh
  sbatch linear_vel_grid.sh

  N_grid=1375
  export N_grid
  sed "s/SIMNAME/${simname}/" linear_vel_grid_template.sh > linear_vel_grid.sh
  sed -i "s/MEMORY/120G/" linear_vel_grid.sh
  sed -i "s/BOXID/${boxid}/" linear_vel_grid.sh
  sed -i "s/N_GRID/${N_grid}/" linear_vel_grid.sh
  sed -i "s/OUTPUT_DIR/test/" linear_vel_grid.sh
  sbatch linear_vel_grid.sh

  done

# 40 GB for 550
# 120 GB for 1375
# Same time for both
