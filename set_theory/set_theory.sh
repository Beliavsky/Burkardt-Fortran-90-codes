#! /bin/bash
#
gfortran -c -Wall set_theory.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv set_theory.o ~/lib/set_theory.o
#
echo "Normal end of execution."
