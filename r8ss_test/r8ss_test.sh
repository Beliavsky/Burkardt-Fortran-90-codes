#! /bin/bash
#
gfortran -c -Wall r8ss_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8ss_test r8ss_test.o $HOME/lib/r8ss.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8ss_test.o
#
./r8ss_test > r8ss_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8ss_test
#
echo "Normal end of execution."
