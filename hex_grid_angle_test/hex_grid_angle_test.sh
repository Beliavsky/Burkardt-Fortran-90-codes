#! /bin/bash
#
gfortran -c -Wall hex_grid_angle_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o hex_grid_angle_test hex_grid_angle_test.o $HOME/lib/hex_grid_angle.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm hex_grid_angle_test.o
#
./hex_grid_angle_test > hex_grid_angle_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm hex_grid_angle_test
#
echo "Normal end of execution."
