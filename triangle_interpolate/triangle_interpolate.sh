#! /bin/bash
#
gfortran -c -Wall triangle_interpolate.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv triangle_interpolate.o ~/lib/triangle_interpolate.o
#
echo "Normal end of execution."
