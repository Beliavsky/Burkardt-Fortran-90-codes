#! /bin/bash
#
gfortran -c -Wall r83s.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r83s.o ~/lib/r83s.o
#
echo "Normal end of execution."
