#! /bin/bash
#
gfortran -c -Wall sgmga.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv sgmga.o ~/lib/sgmga.o
#
echo "Normal end of execution."
