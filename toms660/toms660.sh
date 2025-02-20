#! /bin/bash
#
gfortran -c -Wall toms660.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms660.o ~/lib/toms660.o
#
echo "Normal end of execution."
