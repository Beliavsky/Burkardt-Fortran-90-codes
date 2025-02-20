#! /bin/bash
#
gfortran -c -Wall zero.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv zero.o ~/lib/zero.o
#
echo "Normal end of execution."
