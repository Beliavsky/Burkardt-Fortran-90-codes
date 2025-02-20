#! /bin/bash
#
gfortran -c -Wall burgers_exact_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o burgers_exact_test burgers_exact_test.o $HOME/lib/burgers_exact.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm burgers_exact_test.o
#
./burgers_exact_test > burgers_exact_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm burgers_exact_test
#
echo "Normal end of execution."
