#! /bin/bash
#
gfortran -c -Wall trapezoidal.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv trapezoidal.o ~/lib/trapezoidal.o
#
echo "Normal end of execution."
