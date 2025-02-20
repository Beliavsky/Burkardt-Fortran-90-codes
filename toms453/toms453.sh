#! /bin/bash
#
gfortran -c -Wall toms453.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms453.o ~/lib/toms453.o
#
echo "Normal end of execution."
