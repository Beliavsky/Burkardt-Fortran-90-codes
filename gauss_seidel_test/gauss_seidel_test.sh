#! /bin/bash
#
gfortran -c -Wall gauss_seidel_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran gauss_seidel_test.o $HOME/lib/gauss_seidel.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm gauss_seidel_test.o
#
mv a.out gauss_seidel_test
./gauss_seidel_test > gauss_seidel_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm gauss_seidel_test
#
echo "Normal end of execution."
