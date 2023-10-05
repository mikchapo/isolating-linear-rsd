simname=AbacusCosmos_1100box_planck
fft_smooth=False
export simname
export fft_smooth


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
    sed "s/SIMNAME/${simname}/" vel_grid_smoothing_template.sh > vel_grid_smoothing.sh
    sed -i "s/HOURS/24/" vel_grid_smoothing.sh
    sed -i "s/MEMORY/160G/" vel_grid_smoothing.sh
    sed -i "s/BOXID/${boxid}/" vel_grid_smoothing.sh
    sed -i "s/N_GRID/${N_grid}/" vel_grid_smoothing.sh
    sed -i "s/SMOOTH_TYPE/${smooth_type}/" vel_grid_smoothing.sh
    sed -i "s/R_SMOOTH/${R_smooth}/" vel_grid_smoothing.sh
    sed -i "s/FFT_SMOOTH/${fft_smooth}/" vel_grid_smoothing.sh
    sed -i "s/OUTPUT_DIR/test/" vel_grid_smoothing.sh
    sbatch vel_grid_smoothing.sh
    done

  N_grid=1100
  export N_grid

  for R_smooth in 3 7
    do
    export R_smooth
    sed "s/SIMNAME/${simname}/" vel_grid_smoothing_template.sh > vel_grid_smoothing.sh
    sed -i "s/HOURS/24/" vel_grid_smoothing.sh
    sed -i "s/MEMORY/80G/" vel_grid_smoothing.sh
    sed -i "s/BOXID/${boxid}/" vel_grid_smoothing.sh
    sed -i "s/N_GRID/${N_grid}/" vel_grid_smoothing.sh
    sed -i "s/SMOOTH_TYPE/${smooth_type}/" vel_grid_smoothing.sh
    sed -i "s/R_SMOOTH/${R_smooth}/" vel_grid_smoothing.sh
    sed -i "s/FFT_SMOOTH/${fft_smooth}/" vel_grid_smoothing.sh
    sed -i "s/OUTPUT_DIR/test/" vel_grid_smoothing.sh
    sbatch vel_grid_smoothing.sh
    done

  smooth_type=gaussian
  R_smooth=2
  export smooth_type
  export R_smooth
  sed "s/SIMNAME/${simname}/" vel_grid_smoothing_template.sh > vel_grid_smoothing.sh
  sed -i "s/HOURS/24/" vel_grid_smoothing.sh
  sed -i "s/MEMORY/80G/" vel_grid_smoothing.sh
  sed -i "s/BOXID/${boxid}/" vel_grid_smoothing.sh
  sed -i "s/N_GRID/${N_grid}/" vel_grid_smoothing.sh
  sed -i "s/SMOOTH_TYPE/${smooth_type}/" vel_grid_smoothing.sh
  sed -i "s/R_SMOOTH/${R_smooth}/" vel_grid_smoothing.sh
  sed -i "s/FFT_SMOOTH/${fft_smooth}/" vel_grid_smoothing.sh
  sed -i "s/OUTPUT_DIR/test/" vel_grid_smoothing.sh
  sbatch vel_grid_smoothing.sh
  done
