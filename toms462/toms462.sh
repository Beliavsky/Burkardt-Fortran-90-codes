#! /bin/bash
#
gfortran -c -Wall toms462.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms462.o ~/lib/toms462.o
#
echo "Normal end of execution."
