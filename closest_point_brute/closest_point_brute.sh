#! /bin/bash
#
gfortran -c -Wall closest_point_brute.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv closest_point_brute.o ~/lib/closest_point_brute.o
#
echo "Normal end of execution."
