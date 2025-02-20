#! /bin/bash
#
gfortran -c -Wall r8pbu.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8pbu.o ~/lib/r8pbu.o
#
echo "Normal end of execution."
