#! /bin/bash
#
gfortran -c -Wall r83s_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r83s_test r83s_test.o $HOME/lib/r83s.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r83s_test.o
#
./r83s_test > r83s_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r83s_test
#
echo "Normal end of execution."
