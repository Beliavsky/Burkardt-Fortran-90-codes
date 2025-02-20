#! /bin/bash
#
gfortran -c -Wall knapsack_random.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv knapsack_random.o ~/lib/knapsack_random.o
#
echo "Normal end of execution."
