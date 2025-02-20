#! /bin/bash
#
gfortran -c -Wall prime.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv prime.o ~/lib/prime.o
#
echo "Normal end of execution."
