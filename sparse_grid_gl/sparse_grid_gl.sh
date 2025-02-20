#! /bin/bash
#
gfortran -c -Wall sparse_grid_gl.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sparse_grid_gl.o ~/lib/sparse_grid_gl.o
#
echo "Normal end of execution."
