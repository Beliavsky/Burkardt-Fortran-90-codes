#! /bin/bash
#
gfortran -c -Wall burgers_exact.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv burgers_exact.o ~/lib/burgers_exact.o
#
echo "Normal end of execution."
