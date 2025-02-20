#! /bin/bash
#
gfortran -c -Wall r8ss.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8ss.o ~/lib/r8ss.o
#
echo "Normal end of execution."
