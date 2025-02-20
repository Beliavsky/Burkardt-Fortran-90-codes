#! /bin/bash
#
gfortran -c -Wall polynomial_root_bound_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o polynomial_root_bound_test polynomial_root_bound_test.o \
  $HOME/lib/polynomial_root_bound.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm polynomial_root_bound_test.o
#
./polynomial_root_bound_test > polynomial_root_bound_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm polynomial_root_bound_test
#
echo "Normal end of execution."
