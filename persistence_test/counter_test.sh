#! /bin/bash
#
gfortran -c -Wall counter_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o counter_test counter_test.o /$HOME/lib/counter.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm counter_test.o
#
./counter_test > counter_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm counter_test
#
echo "Normal end of execution."
