#! /bin/bash
#
gfortran -c -Wall r8blt.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8blt.o ~/lib/r8blt.o
#
echo "Normal end of execution."
