#! /bin/bash
#
gfortran -c -Wall r8bto.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8bto.o ~/lib/r8bto.o
#
echo "Normal end of execution."
