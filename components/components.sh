#! /bin/bash
#
gfortran -c -Wall components.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv components.o ~/lib/components.o
#
echo "Normal end of execution."
