#! /bin/bash
#
gfortran -c -Wall r8bb.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8bb.o ~/lib/r8bb.o
#
echo "Normal end of execution."
