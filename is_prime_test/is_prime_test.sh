#! /bin/bash
#
gfortran -c -Wall is_prime_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o is_prime_test is_prime_test.o $HOME/lib/is_prime.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm is_prime_test.o
#
./is_prime_test > is_prime_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm is_prime_test
#
echo "Normal end of execution."
