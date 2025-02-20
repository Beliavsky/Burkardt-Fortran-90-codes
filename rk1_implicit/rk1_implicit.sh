#! /bin/bash
#
gfortran -c -Wall rk1_implicit.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv rk1_implicit.o ~/lib/rk1_implicit.o
#
echo "Normal end of execution."
