#! /bin/bash
#
gfortran -c -Wall toms611.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms611.o ~/lib/toms611.o
#
echo "Normal end of execution."
