#! /bin/bash
#
gfortran -c -Wall prime_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o prime_test prime_test.o $HOME/lib/prime.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm prime_test.o
#
./prime_test > prime_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm prime_test
#
echo "Normal end of execution."
