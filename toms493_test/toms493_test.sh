#! /bin/bash
#
gfortran -c -Wall toms493_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms493_test toms493_test.o $HOME/lib/toms493.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms493_test.o
#
./toms493_test > toms493_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms493_test
#
echo "Normal end of execution."
