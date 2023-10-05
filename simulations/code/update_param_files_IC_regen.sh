for i in 02 13 17 19 20 22 23 24 25 26 29 30 31 32 34 35 36 38
  do
  cd ~/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_${i}_products/
#   sed "s%ic_ngc_nplt_z49.0%ic_ngc_nplt_z49.0_v2%" abacus_ngc_nplt.par > abacus_ngc_nplt_v2.par
  mv ic_ngc_nplt_z49.0 ic_ngc_nplt_z49.0_err
  mv ic_ngc_nplt_z49.0_v2 ic_ngc_nplt_z49.0
  done
