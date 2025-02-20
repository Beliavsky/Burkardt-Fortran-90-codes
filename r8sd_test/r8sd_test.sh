#! /bin/bash
#
gfortran -c -Wall r8sd_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8sd_test r8sd_test.o $HOME/lib/r8sd.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8sd_test.o
#
./r8sd_test > r8sd_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8sd_test
#
echo "Normal end of execution."
