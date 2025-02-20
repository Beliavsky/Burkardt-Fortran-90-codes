#! /bin/bash
#
gfortran -c -Wall r8blt_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8blt_test r8blt_test.o $HOME/lib/r8blt.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8blt_test.o
#
./r8blt_test > r8blt_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8blt_test
#
echo "Normal end of execution."
