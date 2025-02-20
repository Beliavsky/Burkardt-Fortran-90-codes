#! /bin/bash
#
gfortran -c -g -Wall triangle.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv triangle.o ~/lib/triangle.o
#
echo "Normal end of execution."
