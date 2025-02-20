#! /bin/bash
#
gfortran -c -Wall lambert_w.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv lambert_w.o ~/lib/lambert_w.o
#
echo "Normal end of execution."
