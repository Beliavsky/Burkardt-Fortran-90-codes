#! /bin/bash
#
gfortran -c -Wall r8cb.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8cb.o ~/lib/r8cb.o
#
echo "Normal end of execution."
