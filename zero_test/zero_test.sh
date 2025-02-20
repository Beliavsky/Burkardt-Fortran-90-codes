#! /bin/bash
#
gfortran -c -Wall zero_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran zero_test.o $HOME/lib/zero.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm zero_test.o
#
mv a.out zero_test
./zero_test > zero_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm zero_test
#
echo "Normal end of execution."
