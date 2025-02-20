#! /bin/bash
#
gfortran -c -Wall knapsack_brute.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv knapsack_brute.o ~/lib/knapsack_brute.o
#
echo "Normal end of execution."
