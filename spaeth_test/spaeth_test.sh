#! /bin/bash
#
gfortran -c -Wall spaeth_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o spaeth_test spaeth_test.o $HOME/lib/spaeth.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm spaeth_test.o
#
./spaeth_test > spaeth_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm spaeth_test
#
echo "Normal end of execution."
