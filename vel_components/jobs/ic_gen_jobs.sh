simname=AbacusCosmos_1100box_planck
export simname
for boxid in 00-1 00-2 00-3 00-4 00-5 00-6 00-7 00-8 00-9 00-10 00-11 00-12 00-13 00-14 00-15 00-16 00-17 00-18 00-19
# for boxid in 00-0
  do
  export boxid
  sed "s/SIMNAME/${simname}/" ic_gen_template.sh > ic_gen.sh
  sed -i "s/BOXID/${boxid}/" ic_gen.sh
  sbatch ic_gen.sh
  done
