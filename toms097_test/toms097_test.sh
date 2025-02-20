#! /bin/bash
#
gfortran -c -Wall toms097_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms097_test toms097_test.o $HOME/lib/toms097.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms097_test.o
#
./toms097_test > toms097_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms097_test
#
echo "Normal end of execution."
