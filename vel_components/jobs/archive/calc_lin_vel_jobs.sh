simname=AbacusCosmos_1100box_planck
halotype=rockstar
subhalos=True
export simname
export halotype
export subhalos
# for boxid in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
for boxid in 00-0
  do
  export boxid
  sed "s/SIMNAME/${simname}/" calc_lin_vel_template.sh > calc_lin_vel.sh
  sed -i "s/HALOTYPE/${halotype}/" calc_lin_vel.sh
  sed -i "s/BOXID/${boxid}/" calc_lin_vel.sh
  sed -i "s/OUTPUT_DIR/${simname}_${halotype}/" calc_lin_vel.sh
  sed -i "s/SUBHALOS/${subhalos}/" calc_lin_vel.sh
  sbatch calc_lin_vel.sh
  done
