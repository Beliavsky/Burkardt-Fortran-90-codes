#! /bin/bash
#
gfortran -c -Wall middle_square_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o middle_square_test middle_square_test.o $HOME/lib/middle_square.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm middle_square_test.o
#
./middle_square_test > middle_square_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm middle_square_test
#
echo "Normal end of execution."
