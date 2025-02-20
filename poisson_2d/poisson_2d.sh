#! /bin/bash
#
gfortran -c -Wall poisson_2d.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv poisson_2d.o ~/lib/poisson_2d.o
#
echo "Normal end of execution."
