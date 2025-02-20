#! /bin/bash
#
gfortran -c -Wall knapsack_rational.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv knapsack_rational.o ~/lib/knapsack_rational.o
#
echo "Normal end of execution."
