#! /bin/bash
#
gfortran -c -Wall toms446_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms446_test toms446_test.o $HOME/lib/toms446.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms446_test.o
#
./toms446_test > toms446_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms446_test
#
echo "Normal end of execution."
