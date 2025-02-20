#! /bin/bash
#
gfortran -c -Wall r8bb_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8bb_test r8bb_test.o $HOME/lib/r8bb.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8bb_test.o
#
./r8bb_test > r8bb_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8bb_test
#
echo "Normal end of execution."
