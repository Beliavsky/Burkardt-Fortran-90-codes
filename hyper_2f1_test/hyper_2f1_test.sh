#! /bin/bash
#
gfortran -c -Wall hyper_2f1_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o hyper_2f1_test hyper_2f1_test.o \
  $HOME/lib/hyper_2f1.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm hyper_2f1_test.o
#
./hyper_2f1_test > hyper_2f1_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm hyper_2f1_test
#
echo "Normal end of execution."
