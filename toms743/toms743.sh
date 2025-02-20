#! /bin/bash
#
gfortran -c -Wall toms743.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms743.o ~/lib/toms743.o
#
echo "Normal end of execution."
