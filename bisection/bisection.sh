#! /bin/bash
#
gfortran -c -Wall bisection.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv bisection.o ~/lib/bisection.o
#
echo "Normal end of execution."
