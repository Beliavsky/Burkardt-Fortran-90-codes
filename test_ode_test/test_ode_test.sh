#! /bin/bash
#
gfortran -c -Wall test_ode_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o test_ode_test test_ode_test.o $HOME/lib/test_ode.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm test_ode_test.o
#
./test_ode_test > test_ode_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm test_ode_test
#
echo "Normal end of execution."
