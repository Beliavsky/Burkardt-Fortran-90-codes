#! /bin/bash
#
gfortran -c -Wall hex_grid_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o hex_grid_test hex_grid_test.o $HOME/lib/hex_grid.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm hex_grid_test.o
#
./hex_grid_test > hex_grid_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm hex_grid_test
#
echo "Normal end of execution."
