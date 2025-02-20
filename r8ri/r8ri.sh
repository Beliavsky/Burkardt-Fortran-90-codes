#! /bin/bash
#
gfortran -c -Wall r8ri.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8ri.o ~/lib/r8ri.o
#
echo "Normal end of execution."
