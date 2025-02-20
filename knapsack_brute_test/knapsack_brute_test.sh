#! /bin/bash
#
gfortran -c -Wall knapsack_brute_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o knapsack_brute_test knapsack_brute_test.o $HOME/lib/knapsack_brute.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm knapsack_brute_test.o
#
./knapsack_brute_test > knapsack_brute_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm knapsack_brute_test
#
echo "Normal end of execution."
