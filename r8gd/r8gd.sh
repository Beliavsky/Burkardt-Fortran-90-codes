#! /bin/bash
#
gfortran -c -Wall r8gd.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8gd.o ~/lib/r8gd.o
#
echo "Normal end of execution."
