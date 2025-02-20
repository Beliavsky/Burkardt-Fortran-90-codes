#! /bin/bash
#
gfortran -c -Wall levenshtein_distance.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv levenshtein_distance.o ~/lib/levenshtein_distance.o
#
echo "Normal end of execution."
