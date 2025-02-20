#! /bin/bash
#
gfortran -c -Wall r8but_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8but_test r8but_test.o $HOME/lib/r8but.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8but_test.o
#
./r8but_test > r8but_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8but_test
#
echo "Normal end of execution."
