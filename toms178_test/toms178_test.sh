#! /bin/bash
#
gfortran -c -Wall toms178_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms178_test toms178_test.o $HOME/lib/toms178.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms178_test.o
#
./toms178_test > toms178_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms178_test
#
echo "Normal end of execution."
