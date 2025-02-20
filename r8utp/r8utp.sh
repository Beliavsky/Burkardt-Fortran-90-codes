#! /bin/bash
#
gfortran -c -Wall r8utp.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8utp.o ~/lib/r8utp.o
#
echo "Normal end of execution."
