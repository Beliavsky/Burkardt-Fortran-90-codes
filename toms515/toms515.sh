#! /bin/bash
#
gfortran -c -Wall toms515.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms515.o ~/lib/toms515.o
#
echo "Normal end of execution."
