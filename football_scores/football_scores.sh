#! /bin/bash
#
gfortran -c -Wall football_scores.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv football_scores.o ~/lib/football_scores.o
#
echo "Normal end of execution."
