#! /bin/bash
#
gfortran -c -Wall toms757.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms757.o ~/lib/toms757.o
#
echo "Normal end of execution."
