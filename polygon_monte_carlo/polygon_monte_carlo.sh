#! /bin/bash
#
gfortran -c -g -Wall polygon_monte_carlo.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv polygon_monte_carlo.o ~/lib/polygon_monte_carlo.o
#
echo "Normal end of execution."
