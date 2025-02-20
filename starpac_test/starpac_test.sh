#! /bin/bash
#
gfortran -c -fallow-argument-mismatch -Wall starpac_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o starpac_test starpac_test.o $HOME/lib/starpac.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm starpac_test.o
#
./starpac_test > starpac_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm starpac_test
#
echo "Normal end of execution."
