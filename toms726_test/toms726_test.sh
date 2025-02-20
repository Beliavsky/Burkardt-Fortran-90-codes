#! /bin/bash
#
gfortran -c -Wall toms726_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms726_test toms726_test.o $HOME/lib/toms726.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms726_test.o
#
./toms726_test > toms726_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms726_test
#
echo "Normal end of execution."
