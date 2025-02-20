#! /bin/bash
#
gfortran -c -Wall knapsack_values.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv knapsack_values.o ~/lib/knapsack_values.o
#
echo "Normal end of execution."
