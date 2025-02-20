#! /bin/bash
#
gfortran -c -Wall r8pbl_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8pbl_test r8pbl_test.o $HOME/lib/r8pbl.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8pbl_test.o
#
./r8pbl_test > r8pbl_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8pbl_test
#
echo "Normal end of execution."
