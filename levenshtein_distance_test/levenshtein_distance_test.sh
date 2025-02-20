#! /bin/bash
#
gfortran -c -Wall levenshtein_distance_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o levenshtein_distance_test levenshtein_distance_test.o $HOME/lib/levenshtein_distance.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm levenshtein_distance_test.o
#
./levenshtein_distance_test > levenshtein_distance_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm levenshtein_distance_test
#
echo "Normal end of execution."
