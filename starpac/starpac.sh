#! /bin/bash
#
gfortran -c -fallow-argument-mismatch -Wall starpac.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv starpac.o ~/lib/starpac.o
#
echo "Normal end of execution."
