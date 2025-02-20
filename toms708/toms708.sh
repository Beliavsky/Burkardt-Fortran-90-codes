#! /bin/bash
#
gfortran -c -Wall toms708.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms708.o ~/lib/toms708.o
#
echo "Normal end of execution."
