#! /bin/bash
#
gfortran -c -Wall biharmonic_exact.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv biharmonic_exact.o ~/lib/biharmonic_exact.o
#
echo "Normal end of execution."
