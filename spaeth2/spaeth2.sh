#! /bin/bash
#
gfortran -c -Wall spaeth2.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv spaeth2.o ~/lib/spaeth2.o
#
echo "Normal end of execution."
