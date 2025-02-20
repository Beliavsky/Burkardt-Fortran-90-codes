#! /bin/bash
#
gfortran -c -Wall toms515_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms515_test toms515_test.o $HOME/lib/toms515.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms515_test.o
#
./toms515_test > toms515_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms515_test
#
echo "Normal end of execution."
