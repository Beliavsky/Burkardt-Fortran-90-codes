#! /bin/bash
#
gfortran -c -g -Wall rk1.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv rk1.o ~/lib/rk1.o
#
echo "Normal end of execution."
