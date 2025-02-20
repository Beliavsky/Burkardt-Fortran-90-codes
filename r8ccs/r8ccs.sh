#! /bin/bash
#
gfortran -c -Wall r8ccs.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8ccs.o ~/lib/r8ccs.o
#
echo "Normal end of execution."
