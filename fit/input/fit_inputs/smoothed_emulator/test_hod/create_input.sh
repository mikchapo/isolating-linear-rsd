for VOL_FACTOR in 5 20 60
  do
  for RUN_NAME in default rm_small rm_large large_only intermediate_only
    do
    cp test_hod_0_vol-factor-eBOSS_${RUN_NAME}.yaml test_hod_0_vol-factor-${VOL_FACTOR}_${RUN_NAME}.yaml
    sed -i "s/eBOSS/${VOL_FACTOR}/g" test_hod_0_vol-factor-${VOL_FACTOR}_${RUN_NAME}.yaml
    done
  done


for ID in {1..9}
  do
  for VOL_FACTOR in eBOSS 5 20 60
    do
    for RUN_NAME in default rm_small rm_large large_only intermediate_only
      do
      cp test_hod_0_vol-factor-${VOL_FACTOR}_${RUN_NAME}.yaml test_hod_${ID}_vol-factor-${VOL_FACTOR}_${RUN_NAME}.yaml
      sed -i "s/HOD_0/HOD_${ID}/g" test_hod_${ID}_vol-factor-${VOL_FACTOR}_${RUN_NAME}.yaml
      sed -i "s/hod_0/hod_${ID}/g" test_hod_${ID}_vol-factor-${VOL_FACTOR}_${RUN_NAME}.yaml
      done
    done
  done
