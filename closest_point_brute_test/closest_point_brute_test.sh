#! /bin/bash
#
gfortran -c -cpp -Wall closest_point_brute_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o closest_point_brute_test \
  closest_point_brute_test.o \
  $HOME/lib/closest_point_brute.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm closest_point_brute_test.o
#
./closest_point_brute_test > closest_point_brute_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm closest_point_brute_test
#
echo "Normal end of execution."
