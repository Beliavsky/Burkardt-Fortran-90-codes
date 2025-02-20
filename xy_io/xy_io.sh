#! /bin/bash
#
gfortran -c -Wall xy_io.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv xy_io.o ~/lib/xy_io.o
#
echo "Normal end of execution."
