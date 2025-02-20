#! /bin/bash
#
gfortran -c -Wall owen_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o owen_test owen_test.o $HOME/lib/owen.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm owen_test.o
#
./owen_test > owen_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm owen_test
#
echo "Normal end of execution."
