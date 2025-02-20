#! /bin/bash
#
gfortran -c -Wall r8to.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8to.o ~/lib/r8to.o
#
echo "Normal end of execution."
