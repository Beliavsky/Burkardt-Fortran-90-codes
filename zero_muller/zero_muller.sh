#! /bin/bash
#
gfortran -c -Wall zero_muller.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv zero_muller.o ~/lib/zero_muller.o
#
echo "Normal end of execution."
