#! /bin/bash
#
gfortran -c -Wall test_matrix.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv test_matrix.o ~/lib/test_matrix.o
#
echo "Normal end of execution."
