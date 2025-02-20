#! /bin/bash
#
gfortran -c -Wall eigs.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv eigs.o ~/lib/eigs.o
#
echo "Normal end of execution."
