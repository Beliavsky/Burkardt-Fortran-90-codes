#! /bin/bash
#
gfortran -c -Wall hypersphere.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv hypersphere.o ~/lib/hypersphere.o
#
echo "Normal end of execution."
