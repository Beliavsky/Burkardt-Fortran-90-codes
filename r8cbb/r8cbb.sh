#! /bin/bash
#
gfortran -c -Wall r8cbb.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8cbb.o ~/lib/r8cbb.o
#
echo "Normal end of execution."
