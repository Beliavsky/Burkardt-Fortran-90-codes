#! /bin/bash
#
gfortran -c -Wall toms757_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms757_test toms757_test.o $HOME/lib/toms757.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms757_test.o
#
./toms757_test > toms757_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms757_test
#
echo "Normal end of execution."
