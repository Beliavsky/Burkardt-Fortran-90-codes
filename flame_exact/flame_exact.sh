#! /bin/bash
#
gfortran -c -Wall flame_exact.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv flame_exact.o ~/lib/flame_exact.o
#
echo "Normal end of execution."
