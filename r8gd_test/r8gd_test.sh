#! /bin/bash
#
gfortran -c -Wall r8gd_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8gd_test r8gd_test.o $HOME/lib/r8gd.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8gd_test.o
#
./r8gd_test > r8gd_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8gd_test
#
echo "Normal end of execution."
