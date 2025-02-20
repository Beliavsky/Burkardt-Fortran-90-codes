#! /bin/bash
#
gfortran -c -Wall zero_brent_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran zero_brent_test.o $HOME/lib/zero_brent.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm zero_brent_test.o
#
mv a.out zero_brent_test
./zero_brent_test > zero_brent_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm zero_brent_test
#
echo "Normal end of execution."
