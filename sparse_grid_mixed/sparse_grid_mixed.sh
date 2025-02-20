#! /bin/bash
#
gfortran -c -Wall sparse_grid_mixed.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sparse_grid_mixed.o ~/lib/sparse_grid_mixed.o
#
echo "Normal end of execution."
