#! /bin/bash
#
gfortran -c -Wall sparse_grid_hermite.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sparse_grid_hermite.o ~/lib/sparse_grid_hermite.o
#
echo "Normal end of execution."
