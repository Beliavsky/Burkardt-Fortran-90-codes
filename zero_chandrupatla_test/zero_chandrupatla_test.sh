#! /bin/bash
#
gfortran -c -Wall zero_chandrupatla_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o zero_chandrupatla_test zero_chandrupatla_test.o $HOME/lib/zero_chandrupatla.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm zero_chandrupatla_test.o
#
./zero_chandrupatla_test > zero_chandrupatla_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm zero_chandrupatla_test
#
echo "Normal end of execution."
