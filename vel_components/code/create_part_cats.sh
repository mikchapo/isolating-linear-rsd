for BOXID in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
  do
  echo "Starting box ${BOXID}"
  python create_part_cat.py $BOXID 100000
  done
