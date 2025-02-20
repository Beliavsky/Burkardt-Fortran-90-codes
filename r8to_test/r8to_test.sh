#! /bin/bash
#
gfortran -c -Wall r8to_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8to_test r8to_test.o $HOME/lib/r8to.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8to_test.o
#
./r8to_test > r8to_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8to_test
#
echo "Normal end of execution."
