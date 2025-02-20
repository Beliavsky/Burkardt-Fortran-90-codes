#! /bin/bash
#
gfortran -c -Wall poisson_1d_multigrid.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv poisson_1d_multigrid.o ~/lib/poisson_1d_multigrid.o
#
echo "Normal end of execution."
