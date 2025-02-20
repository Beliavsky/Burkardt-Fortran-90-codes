#! /bin/bash
#
gfortran -c -Wall toms243.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms243.o ~/lib/toms243.o
#
echo "Normal end of execution."
