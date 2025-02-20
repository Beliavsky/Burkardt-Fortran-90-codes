#! /bin/bash
#
gfortran -c -Wall gauss_seidel.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv gauss_seidel.o ~/lib/gauss_seidel.o
#
echo "Normal end of execution."
