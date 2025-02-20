#! /bin/bash
#
gfortran -c -Wall normal.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv normal.o ~/lib/normal.o
#
echo "Normal end of execution."
