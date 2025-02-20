#! /bin/bash
#
gfortran -c -Wall toms526.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms526.o ~/lib/toms526.o
#
echo "Normal end of execution."
