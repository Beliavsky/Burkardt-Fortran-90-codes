#! /bin/bash
#
gfortran -c -Wall zero_muller_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o zero_muller_test zero_muller_test.o $HOME/lib/zero_muller.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm zero_muller_test.o
#
./zero_muller_test > zero_muller_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm zero_muller_test
#
echo "Normal end of execution."
