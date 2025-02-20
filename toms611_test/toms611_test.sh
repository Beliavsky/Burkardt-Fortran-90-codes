#! /bin/bash
#
gfortran -c -Wall toms611_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms611_test toms611_test.o $HOME/lib/toms611.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms611_test.o
#
./toms611_test > toms611_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms611_test
#
echo "Normal end of execution."
