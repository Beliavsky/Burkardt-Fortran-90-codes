#! /bin/bash
#
gfortran -c -Wall -fopenmp helmholtz_openmp.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -fopenmp helmholtz_openmp.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm helmholtz_openmp.o
mv a.out $HOME/bin/helmholtz_openmp
#
echo "Normal end of execution."
