#! /bin/bash
#
gfortran -c -Wall zero_itp.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv zero_itp.o ~/lib/zero_itp.o
#
echo "Normal end of execution."
