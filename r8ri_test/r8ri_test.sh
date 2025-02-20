#! /bin/bash
#
gfortran -c -Wall r8ri_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8ri_test r8ri_test.o $HOME/lib/r8ri.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8ri_test.o
#
./r8ri_test > r8ri_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8ri_test
#
echo "Normal end of execution."
