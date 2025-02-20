#! /bin/bash
#
gfortran -c -Wall r8sm.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8sm.o ~/lib/r8sm.o
#
echo "Normal end of execution."
