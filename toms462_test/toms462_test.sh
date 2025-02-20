#! /bin/bash
#
gfortran -c -Wall toms462_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms462_test toms462_test.o $HOME/lib/toms462.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms462_test.o
#
./toms462_test > toms462_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms462_test
#
echo "Normal end of execution."
