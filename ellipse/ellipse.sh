#! /bin/bash
#
gfortran -c -Wall ellipse.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv ellipse.o ~/lib/ellipse.o
#
echo "Normal end of execution."
