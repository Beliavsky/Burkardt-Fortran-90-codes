#! /bin/bash
#
gfortran -c -Wall subset_sum_swap.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv subset_sum_swap.o ~/lib/subset_sum_swap.o
#
echo "Normal end of execution."
