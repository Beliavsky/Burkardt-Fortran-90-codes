#! /bin/bash
#
gfortran -c -Wall toms097.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms097.o ~/lib/toms097.o
#
echo "Normal end of execution."
