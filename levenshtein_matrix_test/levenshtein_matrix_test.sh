#! /bin/bash
#
gfortran -c -Wall levenshtein_matrix_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o levenshtein_matrix_test levenshtein_matrix_test.o $HOME/lib/levenshtein_matrix.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm levenshtein_matrix_test.o
#
./levenshtein_matrix_test > levenshtein_matrix_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm levenshtein_matrix_test
#
echo "Normal end of execution."
