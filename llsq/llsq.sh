#! /bin/bash
#
gfortran -c -Wall llsq.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv llsq.o ~/lib/llsq.o
#
echo "Normal end of execution."
