#! /bin/bash
#
gfortran -c -Wall zero_laguerre.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv zero_laguerre.o ~/lib/zero_laguerre.o
#
echo "Normal end of execution."
