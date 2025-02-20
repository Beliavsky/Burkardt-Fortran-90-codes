#! /bin/bash
#
gfortran -c -Wall toms661_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms661_test toms661_test.o $HOME/lib/toms661.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms661_test.o
#
./toms661_test > toms661_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms661_test
#
echo "Normal end of execution."
