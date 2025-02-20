#! /bin/bash
#
gfortran -c -Wall r83p.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r83p.o ~/lib/r83p.o
#
echo "Normal end of execution."
