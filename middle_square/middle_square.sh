#! /bin/bash
#
gfortran -c -Wall middle_square.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv middle_square.o ~/lib/middle_square.o
#
echo "Normal end of execution."
