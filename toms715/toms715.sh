#! /bin/bash
#
gfortran -c -Wall toms715.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms715.o ~/lib/toms715.o
#
echo "Normal end of execution."
