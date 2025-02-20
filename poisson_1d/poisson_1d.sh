#! /bin/bash
#
gfortran -c -Wall poisson_1d.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv poisson_1d.o ~/lib/poisson_1d.o
#
echo "Normal end of execution."
