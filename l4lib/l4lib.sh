#! /bin/bash
#
gfortran -c -g -Wall l4lib.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv l4lib.o ~/lib/l4lib.o
#
echo "Normal end of execution."
