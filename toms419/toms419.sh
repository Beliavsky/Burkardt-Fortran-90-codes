#! /bin/bash
#
gfortran -c -Wall toms419.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms419.o ~/lib/toms419.o
#
echo "Normal end of execution."
