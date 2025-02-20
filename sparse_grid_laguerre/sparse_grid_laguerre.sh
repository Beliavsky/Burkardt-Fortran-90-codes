#! /bin/bash
#
gfortran -c -Wall sparse_grid_laguerre.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sparse_grid_laguerre.o ~/lib/sparse_grid_laguerre.o
#
echo "Normal end of execution."
