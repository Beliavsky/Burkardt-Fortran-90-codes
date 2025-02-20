#! /bin/bash
#
gfortran -c -Wall subset_sum_swap_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o subset_sum_swap_test \
  subset_sum_swap_test.o \
  $HOME/lib/subset_sum_swap.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm subset_sum_swap_test.o
#
./subset_sum_swap_test > subset_sum_swap_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm subset_sum_swap_test
#
echo "Normal end of execution."
