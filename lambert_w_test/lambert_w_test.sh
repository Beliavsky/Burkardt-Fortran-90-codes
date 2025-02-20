#! /bin/bash
#
gfortran -c -Wall lambert_w_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o lambert_w_test lambert_w_test.o $HOME/lib/lambert_w.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm lambert_w_test.o
#
./lambert_w_test > lambert_w_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm lambert_w_test
#
echo "Normal end of execution."
