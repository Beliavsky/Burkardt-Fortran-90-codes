#! /bin/bash
#
gfortran -c -Wall hex_grid.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv hex_grid.o ~/lib/hex_grid.o
#
echo "Normal end of execution."
