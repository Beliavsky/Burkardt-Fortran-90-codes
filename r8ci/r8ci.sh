#! /bin/bash
#
gfortran -c -Wall r8ci.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8ci.o ~/lib/r8ci.o
#
echo "Normal end of execution."
