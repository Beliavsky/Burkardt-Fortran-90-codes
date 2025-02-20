#! /bin/bash
#
gfortran -c -Wall midpoint.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv midpoint.o ~/lib/midpoint.o
#
echo "Normal end of execution."
