#! /bin/bash
#
gfortran -c -Wall slap.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv slap.o ~/lib/slap.o
#
echo "Normal end of execution."
