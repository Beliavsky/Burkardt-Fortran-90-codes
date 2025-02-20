#! /bin/bash
#
gfortran -c -Wall toms655_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms655_test toms655_test.o $HOME/lib/toms655.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms655_test.o
#
./toms655_test > toms655_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms655_test
#
echo "Normal end of execution."
