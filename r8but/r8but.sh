#! /bin/bash
#
gfortran -c -Wall r8but.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8but.o ~/lib/r8but.o
#
echo "Normal end of execution."
