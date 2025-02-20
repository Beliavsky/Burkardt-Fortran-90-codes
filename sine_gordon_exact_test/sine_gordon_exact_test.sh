#! /bin/bash
#
gfortran -c -Wall sine_gordon_exact_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o sine_gordon_exact_test sine_gordon_exact_test.o $HOME/lib/sine_gordon_exact.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sine_gordon_exact_test.o
#
./sine_gordon_exact_test > sine_gordon_exact_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sine_gordon_exact_test
#
echo "Normal end of execution."
