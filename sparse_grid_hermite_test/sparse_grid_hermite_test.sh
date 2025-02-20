#! /bin/bash
#
gfortran -c -Wall sparse_grid_hermite_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran sparse_grid_hermite_test.o $HOME/lib/sparse_grid_hermite.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sparse_grid_hermite_test.o
#
mv a.out sparse_grid_hermite_test
./sparse_grid_hermite_test > sparse_grid_hermite_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sparse_grid_hermite_test
#
echo "Normal end of execution."
