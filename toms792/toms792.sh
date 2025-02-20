#! /bin/bash
#
gfortran -c -Wall toms792.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv toms792.o ~/lib/toms792.o
#
echo "Normal end of execution."
