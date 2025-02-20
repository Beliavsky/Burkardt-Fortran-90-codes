#! /bin/bash
#
gfortran -c -Wall bisection_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o bisection_test bisection_test.o $HOME/lib/bisection.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm bisection_test.o
#
./bisection_test > bisection_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm bisection_test
#
echo "Normal end of execution."
