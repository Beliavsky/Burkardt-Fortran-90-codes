#! /bin/bash
#
gfortran -c -Wall r8pbl.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8pbl.o ~/lib/r8pbl.o
#
echo "Normal end of execution."
