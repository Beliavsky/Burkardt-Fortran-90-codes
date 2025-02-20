#! /bin/bash
#
gfortran -c -Wall polynomial_root_bound.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv polynomial_root_bound.o ~/lib/polynomial_root_bound.o
#
echo "Normal end of execution."
