#! /bin/bash
#
gfortran -c -Wall r83t.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv r83t.o ~/lib/r83t.o
#
echo "Normal end of execution."
