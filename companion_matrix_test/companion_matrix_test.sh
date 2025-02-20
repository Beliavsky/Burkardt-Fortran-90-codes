#! /bin/bash
#
gfortran -c -Wall companion_matrix_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o companion_matrix_test \
  companion_matrix_test.o \
  $HOME/lib/companion_matrix.o \
  $HOME/lib/eigs.o \
  $HOME/lib/toms493.o \
  -llapack
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm companion_matrix_test.o
#
./companion_matrix_test > companion_matrix_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm companion_matrix_test
#
echo "Normal end of execution."
