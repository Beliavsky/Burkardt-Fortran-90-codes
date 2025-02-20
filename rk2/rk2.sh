#! /bin/bash
#
gfortran -c -Wall rk2.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv rk2.o ~/lib/rk2.o
#
echo "Normal end of execution."
