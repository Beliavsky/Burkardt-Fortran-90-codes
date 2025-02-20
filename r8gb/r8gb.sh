#! /bin/bash
#
gfortran -c -Wall r8gb.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8gb.o ~/lib/r8gb.o
#
echo "Normal end of execution."
