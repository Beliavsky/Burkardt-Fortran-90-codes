#! /bin/bash
#
gfortran -c -g -Wall unicycle.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv unicycle.o ~/lib/unicycle.o
#
echo "Normal end of execution."
