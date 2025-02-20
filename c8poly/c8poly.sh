#! /bin/bash
#
gfortran -c -Wall c8poly.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv c8poly.o ~/lib
#
echo "Normal end of execution."
