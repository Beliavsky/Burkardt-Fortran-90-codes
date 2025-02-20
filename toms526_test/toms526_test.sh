#! /bin/bash
#
gfortran -c -Wall toms526_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms526_test toms526_test.o $HOME/lib/toms526.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms526_test.o
#
./toms526_test > toms526_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms526_test
#
echo "Normal end of execution."
