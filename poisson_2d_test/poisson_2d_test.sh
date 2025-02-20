#! /bin/bash
#
gfortran -c -Wall poisson_2d_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o poisson_2d_test poisson_2d_test.o $HOME/lib/poisson_2d.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm poisson_2d_test.o
#
./poisson_2d_test > poisson_2d_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm poisson_2d_test
#
echo "Normal end of execution."
