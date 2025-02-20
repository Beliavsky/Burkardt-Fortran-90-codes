#! /bin/bash
#
gfortran -c -Wall r8ncf.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8ncf.o ~/lib/r8ncf.o
#
echo "Normal end of execution."
