#! /bin/bash
#
gfortran -c -Wall sine_gordon_exact.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sine_gordon_exact.o ~/lib/sine_gordon_exact.o
#
echo "Normal end of execution."
