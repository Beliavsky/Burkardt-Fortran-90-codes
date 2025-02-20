#! /bin/bash
#
gfortran -c -Wall steam.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv steam.o ~/lib/steam.o
#
echo "Normal end of execution."
