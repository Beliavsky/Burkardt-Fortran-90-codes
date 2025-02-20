#! /bin/bash
#
gfortran -c -Wall zero_laguerre_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o zero_laguerre_test zero_laguerre_test.o $HOME/lib/zero_laguerre.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm zero_laguerre_test.o
#
./zero_laguerre_test > zero_laguerre_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm zero_laguerre_test
#
echo "Normal end of execution."
