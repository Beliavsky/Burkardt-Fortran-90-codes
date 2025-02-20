#! /bin/bash
#
gfortran -c -Wall zero_brent.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv zero_brent.o ~/lib/zero_brent.o
#
echo "Normal end of execution."
