#! /bin/bash
#
gfortran -c -Wall toms715_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms715_test toms715_test.o $HOME/lib/toms715.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms715_test.o
#
./toms715_test > toms715_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms715_test
#
echo "Normal end of execution."
