#! /bin/bash
#
gfortran -c -Wall toms112.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms112.o ~/lib/toms112.o
#
echo "Normal end of execution."
