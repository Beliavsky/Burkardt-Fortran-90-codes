#! /bin/bash
#
gfortran -c -Wall r8ltt.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r8ltt.o ~/lib/r8ltt.o
#
echo "Normal end of execution."
