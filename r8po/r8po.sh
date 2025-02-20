#! /bin/bash
#
gfortran -c -Wall r8po.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8po.o ~/lib/r8po.o
#
echo "Normal end of execution."
