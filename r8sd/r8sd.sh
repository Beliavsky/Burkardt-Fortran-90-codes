#! /bin/bash
#
gfortran -c -Wall r8sd.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8sd.o ~/lib/r8sd.o
#
echo "Normal end of execution."
