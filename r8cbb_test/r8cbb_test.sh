#! /bin/bash
#
gfortran -c -Wall r8cbb_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8cbb_test r8cbb_test.o $HOME/lib/r8cbb.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8cbb_test.o
#
./r8cbb_test > r8cbb_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8cbb_test
#
echo "Normal end of execution."
