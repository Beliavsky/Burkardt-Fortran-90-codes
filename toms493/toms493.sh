#! /bin/bash
#
gfortran -c -Wall toms493.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms493.o ~/lib/toms493.o
#
echo "Normal end of execution."
