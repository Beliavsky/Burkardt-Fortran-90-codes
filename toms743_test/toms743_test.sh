#! /bin/bash
#
gfortran -c -Wall toms743_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms743_test toms743_test.o $HOME/lib/toms743.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms743_test.o
#
./toms743_test > toms743_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms743_test
#
echo "Normal end of execution."
