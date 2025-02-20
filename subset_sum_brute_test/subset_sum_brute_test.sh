#! /bin/bash
#
gfortran -c -cpp -Wall subset_sum_brute_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o subset_sum_brute_test subset_sum_brute_test.o $HOME/lib/subset_sum_brute.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm subset_sum_brute_test.o
#
./subset_sum_brute_test > subset_sum_brute_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm subset_sum_brute_test
#
echo "Normal end of execution."
