#! /bin/bash
#
gfortran -c -Wall glomin.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv glomin.o ~/lib/glomin.o
#
echo "Normal end of execution."
