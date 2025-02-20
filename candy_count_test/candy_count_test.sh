#! /bin/bash
#
gfortran -c -Wall candy_count_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o candy_count_test candy_count_test.o $HOME/lib/candy_count.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm candy_count_test.o
#
./candy_count_test > candy_count_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm candy_count_test
#
echo "Normal end of execution."
