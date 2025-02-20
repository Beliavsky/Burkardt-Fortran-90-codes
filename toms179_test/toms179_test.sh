#! /bin/bash
#
gfortran -c -Wall toms179_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms179_test toms179_test.o $HOME/lib/toms179.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms179_test.o
#
./toms179_test > toms179_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms179_test
#
echo "Normal end of execution."
