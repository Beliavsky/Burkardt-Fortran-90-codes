#! /bin/bash
#
gfortran -c -Wall r85_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r85_test r85_test.o $HOME/lib/r85.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r85_test.o
#
./r85_test > r85_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r85_test
#
echo "Normal end of execution."
