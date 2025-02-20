#! /bin/bash
#
gfortran -c -Wall knapsack_values_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o knapsack_values_test knapsack_values_test.o $HOME/lib/knapsack_values.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm knapsack_values_test.o
#
./knapsack_values_test > knapsack_values_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm knapsack_values_test
#
echo "Normal end of execution."
