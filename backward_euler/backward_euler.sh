#! /bin/bash
#
gfortran -c -Wall backward_euler.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv backward_euler.o ~/lib/backward_euler.o
#
echo "Normal end of execution."
