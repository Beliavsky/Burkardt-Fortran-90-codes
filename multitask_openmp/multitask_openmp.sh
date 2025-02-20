#! /bin/bash
#
gfortran -c -Wall -fopenmp multitask_openmp.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -fopenmp multitask_openmp.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm multitask_openmp.o
mv a.out $HOME/bin/multitask_openmp
#
echo "Normal end of execution."
