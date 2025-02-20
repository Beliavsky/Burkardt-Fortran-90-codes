#! /bin/bash
#
gfortran -c -Wall xyz_io.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv xyz_io.o ~/lib/xyz_io.o
#
echo "Normal end of execution."
