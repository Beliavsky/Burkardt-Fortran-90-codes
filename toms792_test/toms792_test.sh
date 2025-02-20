#! /bin/bash
#
gfortran -c -Wall toms792_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o toms792_test toms792_test.o $HOME/lib/toms792.o $HOME/lib/toms790.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm toms792_test.o
#
./toms792_test < input.txt > toms792_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm toms792_test
#
echo "Normal end of execution."
