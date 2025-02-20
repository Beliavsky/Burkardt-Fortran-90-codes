#! /bin/bash
#
gfortran -c -Wall rk2_implicit.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv rk2_implicit.o ~/lib/rk2_implicit.o
#
echo "Normal end of execution."
