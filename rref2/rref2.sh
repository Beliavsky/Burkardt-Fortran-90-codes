#! /bin/bash
#
gfortran -c -g -Wall rref2.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv rref2.o ~/lib/rref2.o
#
echo "Normal end of execution."
