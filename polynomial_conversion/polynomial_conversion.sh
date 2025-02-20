#! /bin/bash
#
gfortran -c -Wall polynomial_conversion.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv polynomial_conversion.o ~/lib/polynomial_conversion.o
#
echo "Normal end of execution."
