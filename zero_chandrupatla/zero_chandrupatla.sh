#! /bin/bash
#
gfortran -c -Wall zero_chandrupatla.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv zero_chandrupatla.o ~/lib/zero_chandrupatla.o
#
echo "Normal end of execution."
