#! /bin/bash
#
gfortran -c -Wall r83v.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r83v.o ~/lib/r83v.o
#
echo "Normal end of execution."
