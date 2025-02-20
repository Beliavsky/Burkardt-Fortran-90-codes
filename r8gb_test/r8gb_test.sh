#! /bin/bash
#
gfortran -c -Wall r8gb_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8gb_test r8gb_test.o $HOME/lib/r8gb.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8gb_test.o
#
./r8gb_test > r8gb_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8gb_test
#
echo "Normal end of execution."
