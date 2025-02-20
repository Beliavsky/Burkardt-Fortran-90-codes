#! /bin/bash
#
gfortran -c -Wall toms659.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms659.o ~/lib/toms659.o
#
echo "Normal end of execution."
