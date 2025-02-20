#! /bin/bash
#
gfortran -c -Wall sphere_triangle_monte_carlo.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sphere_triangle_monte_carlo.o ~/lib/sphere_triangle_monte_carlo.o
#
echo "Normal end of execution."
