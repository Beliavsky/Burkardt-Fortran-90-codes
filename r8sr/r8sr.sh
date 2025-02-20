#! /bin/bash
#
gfortran -c -Wall r8sr.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8sr.o ~/lib/r8sr.o
#
echo "Normal end of execution."
