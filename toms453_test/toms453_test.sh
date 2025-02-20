#! /bin/bash
#
gfortran -c -Wall toms453_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms453_test toms453_test.o $HOME/lib/toms453.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms453_test.o
#
./toms453_test > toms453_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms453_test
#
echo "Normal end of execution."
