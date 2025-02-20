#! /bin/bash
#
gfortran -c -Wall spaeth.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv spaeth.o ~/lib/spaeth.o
#
echo "Normal end of execution."
