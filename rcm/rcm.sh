#! /bin/bash
#
gfortran -c -Wall rcm.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv rcm.o ~/lib/rcm.o
#
echo "Normal end of execution."
