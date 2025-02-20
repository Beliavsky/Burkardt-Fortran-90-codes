#! /bin/bash
#
gfortran -c -Wall hypergeometric_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o hypergeometric_test hypergeometric_test.o $HOME/lib/hypergeometric.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm hypergeometric_test.o
#
./hypergeometric_test > hypergeometric_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm hypergeometric_test
#
echo "Normal end of execution."
