#! /bin/bash
#
gfortran -c -Wall closest_pair_brute.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv closest_pair_brute.o ~/lib/closest_pair_brute.o
#
echo "Normal end of execution."
