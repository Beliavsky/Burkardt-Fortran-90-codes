#! /bin/bash
#
gfortran -c -Wall eispack.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv eispack.o ~/lib/eispack.o
#
echo "Normal end of execution."
