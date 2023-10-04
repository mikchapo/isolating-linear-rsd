for VOL_FACTOR in eBOSS 5 60
  do
  for RUN_NAME in default scales
    do
    cp test_hod_0_vol-factor-20_${RUN_NAME}.yaml test_hod_0_vol-factor-${VOL_FACTOR}_${RUN_NAME}.yaml
    sed -i "s/20/${VOL_FACTOR}/g" test_hod_0_vol-factor-${VOL_FACTOR}_${RUN_NAME}.yaml
    done
  done


for ID in {1..9}
  do
  for RUN_NAME in default scales
    do
    cp test_hod_0_vol-factor-20_${RUN_NAME}.yaml test_hod_${ID}_vol-factor-20_${RUN_NAME}.yaml
    sed -i "s/test_hod: 0/test_hod: ${ID}/g" test_hod_${ID}_vol-factor-20_${RUN_NAME}.yaml
    sed -i "s/hod_0/hod_${ID}/g" test_hod_${ID}_vol-factor-20_${RUN_NAME}.yaml
    done
  cp test_hod_0_vol-factors.yaml test_hod_${ID}_vol-factors.yaml
  sed -i "s/test_hod: 0/test_hod: ${ID}/g" test_hod_${ID}_vol-factors.yaml
  sed -i "s/hod_0/hod_${ID}/g" test_hod_${ID}_vol-factors.yaml
  done
