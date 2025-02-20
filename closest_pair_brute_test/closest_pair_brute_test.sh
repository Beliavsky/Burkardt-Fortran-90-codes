#! /bin/bash
#
gfortran -c -Wall closest_pair_brute_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o closest_pair_brute_test closest_pair_brute_test.o $HOME/lib/closest_pair_brute.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm closest_pair_brute_test.o
#
./closest_pair_brute_test > closest_pair_brute_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm closest_pair_brute_test
#
echo "Normal end of execution."
