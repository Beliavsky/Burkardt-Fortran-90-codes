#! /bin/bash
#
gfortran -c -Wall toms886.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms886.o ~/lib/toms886.o
#
echo "Normal end of execution."
