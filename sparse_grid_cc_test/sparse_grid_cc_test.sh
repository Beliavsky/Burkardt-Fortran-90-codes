#! /bin/bash
#
gfortran -c -Wall sparse_grid_cc_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran sparse_grid_cc_test.o $HOME/lib/sparse_grid_cc.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sparse_grid_cc_test.o
#
mv a.out sparse_grid_cc_test
./sparse_grid_cc_test > sparse_grid_cc_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sparse_grid_cc_test
#
echo "Normal end of execution."
