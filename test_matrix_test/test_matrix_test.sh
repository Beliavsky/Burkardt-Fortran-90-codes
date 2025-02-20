#! /bin/bash
#
gfortran -c -Wall test_matrix_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o test_matrix_test test_matrix_test.o /$HOME/lib/test_matrix.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm test_matrix_test.o
#
./test_matrix_test > test_matrix_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm test_matrix_test
#
echo "Normal end of execution."
