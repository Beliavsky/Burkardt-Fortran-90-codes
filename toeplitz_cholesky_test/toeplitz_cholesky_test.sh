#! /bin/bash
#
gfortran -c -Wall toeplitz_cholesky_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran toeplitz_cholesky_test.o $HOME/lib/toeplitz_cholesky.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toeplitz_cholesky_test.o
#
mv a.out toeplitz_cholesky_test
./toeplitz_cholesky_test > toeplitz_cholesky_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toeplitz_cholesky_test
#
echo "Normal end of execution."
