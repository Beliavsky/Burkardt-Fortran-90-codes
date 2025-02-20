#! /bin/bash
#
gfortran -c -Wall r85.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r85.o ~/lib/r85.o
#
echo "Normal end of execution."
