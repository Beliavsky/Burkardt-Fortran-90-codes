#! /bin/bash
#
gfortran -c -Wall hypergeometric.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv hypergeometric.o ~/lib/hypergeometric.o
#
echo "Normal end of execution."
