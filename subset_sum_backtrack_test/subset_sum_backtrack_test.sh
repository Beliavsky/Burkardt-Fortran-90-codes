#! /bin/bash
#
gfortran -c -Wall subset_sum_backtrack_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran subset_sum_backtrack_test.o $HOME/lib/subset_sum_backtrack.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm subset_sum_backtrack_test.o
#
mv a.out subset_sum_backtrack_test
./subset_sum_backtrack_test > subset_sum_backtrack_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm subset_sum_backtrack_test
#
echo "Normal end of execution."
