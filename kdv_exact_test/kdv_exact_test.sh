#! /bin/bash
#
gfortran -c -Wall kdv_exact_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o kdv_exact_test kdv_exact_test.o $HOME/lib/kdv_exact.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm kdv_exact_test.o
#
./kdv_exact_test > kdv_exact_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm kdv_exact_test
#
echo "Normal end of execution."
