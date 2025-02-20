#! /bin/bash
#
gfortran -c -Wall sparse_grid_cc.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sparse_grid_cc.o ~/lib/sparse_grid_cc.o
#
echo "Normal end of execution."
