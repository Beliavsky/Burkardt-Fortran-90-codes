#! /bin/bash
#
gfortran -c -g -Wall rref2_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o rref2_test rref2_test.o $HOME/lib/rref2.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm rref2_test.o
#
./rref2_test > rref2_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm rref2_test
#
echo "Normal end of execution."
