# simname=AbacusCosmos_1100box_planck
simname=AbacusCosmos_1100box
export simname
# for boxid in 00-0 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
for boxid in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 planck
  do
  export boxid
  sed "s/SIMNAME/${simname}/" avoid_purge_template.sh > avoid_purge.sh
  sed -i "s/BOXID/${boxid}/" avoid_purge.sh
  sbatch avoid_purge.sh
  sleep 1
  done
