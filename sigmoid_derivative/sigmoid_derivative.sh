#! /bin/bash
#
gfortran -c -Wall sigmoid_derivative.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sigmoid_derivative.o ~/lib/sigmoid_derivative.o
#
echo "Normal end of execution."
