#! /bin/bash
#
gfortran -c -Wall r8sto.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8sto.o ~/lib/r8sto.o
#
echo "Normal end of execution."
