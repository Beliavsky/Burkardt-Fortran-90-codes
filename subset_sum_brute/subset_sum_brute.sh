#! /bin/bash
#
gfortran -c -Wall subset_sum_brute.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv subset_sum_brute.o ~/lib/subset_sum_brute.o
#
echo "Normal end of execution."
