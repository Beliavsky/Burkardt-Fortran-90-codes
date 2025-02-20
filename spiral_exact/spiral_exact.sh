#! /bin/bash
#
gfortran -c -Wall spiral_exact.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv spiral_exact.o ~/lib/spiral_exact.o
#
echo "Normal end of execution."
