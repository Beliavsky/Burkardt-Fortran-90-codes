#! /bin/bash
#
gfortran -c -Wall r8vm.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8vm.o ~/lib/r8vm.o
#
echo "Normal end of execution."
