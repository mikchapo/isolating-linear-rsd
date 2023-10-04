# for input in eboss_fiducial.yaml eboss_rm_small.yaml eboss_rm_intermediate.yaml eboss_rm_large.yaml eboss_large_only.yaml eboss_mono+quad.yaml eboss_mono+wp.yaml eboss_fixed_gamma_n.yaml eboss_fixed_gamma_l.yaml eboss_equal_gammas.yaml
for input in eboss_rm_small.yaml eboss_rm_intermediate.yaml eboss_rm_large.yaml eboss_large_only.yaml eboss_mono+quad.yaml eboss_mono+wp.yaml
  do
  sed -i "s%chains%chains/smoothed_emulator%" $input
  sed -i "s/abacus_cov_mat/abacus_smooth-vel_cov_mat/" $input
  sed -i "s/aemulus_only/emulator_only/" $input
  done
