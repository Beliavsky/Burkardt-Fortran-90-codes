#! /bin/bash
#
gfortran -c -Wall hex_grid_angle.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv hex_grid_angle.o ~/lib/hex_grid_angle.o
#
echo "Normal end of execution."
