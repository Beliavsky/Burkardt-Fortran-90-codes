#! /bin/bash
#
gfortran -c -g -Wall l4lib_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o l4lib_test l4lib_test.o /$HOME/lib/l4lib.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm l4lib_test.o
#
./l4lib_test > l4lib_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm l4lib_test
#
echo "Normal end of execution."
