#! /bin/bash
#
gfortran -c -Wall quadpack.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv quadpack.o ~/lib/quadpack.o
#
echo "Normal end of execution."
