#! /bin/bash
#
gfortran -c -Wall sparse_count_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o sparse_count_test sparse_count_test.o $HOME/lib/sparse_count.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sparse_count_test.o
#
./sparse_count_test > sparse_count_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sparse_count_test
#
echo "Normal end of execution."
