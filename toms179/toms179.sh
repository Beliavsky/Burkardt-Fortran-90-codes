#! /bin/bash
#
gfortran -c -Wall toms179.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms179.o ~/lib/toms179.o
#
echo "Normal end of execution."
