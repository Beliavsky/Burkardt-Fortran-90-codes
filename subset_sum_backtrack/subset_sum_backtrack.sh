#! /bin/bash
#
gfortran -c -Wall subset_sum_backtrack.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv subset_sum_backtrack.o ~/lib/subset_sum_backtrack.o
#
echo "Normal end of execution."
