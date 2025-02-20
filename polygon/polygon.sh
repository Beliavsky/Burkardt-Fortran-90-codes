#! /bin/bash
#
gfortran -c -g -Wall polygon.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv polygon.o ~/lib/polygon.o
#
echo "Normal end of execution."
