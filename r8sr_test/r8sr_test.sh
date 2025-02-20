#! /bin/bash
#
gfortran -c -Wall r8sr_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8sr_test r8sr_test.o $HOME/lib/r8sr.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8sr_test.o
#
./r8sr_test > r8sr_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8sr_test
#
echo "Normal end of execution."
