#! /bin/bash
#
gfortran -c -Wall levenshtein_matrix.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv levenshtein_matrix.o ~/lib/levenshtein_matrix.o
#
echo "Normal end of execution."
