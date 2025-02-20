#! /bin/bash
#
~/bin/nint_exactness_mixed sparse_grid_mixed_d2_l2_ccxcc 7 > nint_exactness_mixed_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
