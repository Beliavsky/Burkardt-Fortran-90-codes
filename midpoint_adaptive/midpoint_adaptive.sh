#! /bin/bash
#
gfortran -c -Wall midpoint_adaptive.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv midpoint_adaptive.o ~/lib/midpoint_adaptive.o
#
echo "Normal end of execution."
