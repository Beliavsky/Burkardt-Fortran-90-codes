#! /bin/bash
#
gfortran -c -Wall toms112_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms112_test toms112_test.o $HOME/lib/toms112.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms112_test.o
#
./toms112_test > toms112_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms112_test
#
echo "Normal end of execution."
