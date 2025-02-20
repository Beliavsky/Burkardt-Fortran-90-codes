#! /bin/bash
#
gfortran -c -Wall toms655.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms655.o ~/lib/toms655.o
#
echo "Normal end of execution."
