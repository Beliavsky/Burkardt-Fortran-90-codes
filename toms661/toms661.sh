#! /bin/bash
#
gfortran -c -Wall toms661.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms661.o ~/lib/toms661.o
#
echo "Normal end of execution."
