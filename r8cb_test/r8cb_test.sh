#! /bin/bash
#
gfortran -c -Wall r8cb_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8cb_test r8cb_test.o $HOME/lib/r8cb.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8cb_test.o
#
./r8cb_test > r8cb_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8cb_test
#
echo "Normal end of execution."
