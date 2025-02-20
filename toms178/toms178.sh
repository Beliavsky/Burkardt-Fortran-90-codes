#! /bin/bash
#
gfortran -c -Wall toms178.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms178.o ~/lib/toms178.o
#
echo "Normal end of execution."
