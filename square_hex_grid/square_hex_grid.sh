#! /bin/bash
#
gfortran -c -Wall square_hex_grid.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv square_hex_grid.o ~/lib/square_hex_grid.o
#
echo "Normal end of execution."
