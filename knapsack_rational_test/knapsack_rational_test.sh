#! /bin/bash
#
gfortran -c -Wall knapsack_rational_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o knapsack_rational_test knapsack_rational_test.o $HOME/lib/knapsack_rational.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm knapsack_rational_test.o
#
./knapsack_rational_test > knapsack_rational_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm knapsack_rational_test
#
echo "Normal end of execution."
