#! /bin/bash
#
gfortran -c -Wall ellipsoid.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv ellipsoid.o ~/lib/ellipsoid.o
#
echo "Normal end of execution."
