#! /bin/bash
#
gfortran -c -Wall toms446.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms446.o ~/lib/toms446.o
#
echo "Normal end of execution."
