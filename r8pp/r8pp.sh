#! /bin/bash
#
gfortran -c -Wall r8pp.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8pp.o ~/lib/r8pp.o
#
echo "Normal end of execution."
