#! /bin/bash
#
gfortran -c -Wall poisson_1d_multigrid_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran poisson_1d_multigrid_test.o $HOME/lib/poisson_1d_multigrid.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm poisson_1d_multigrid_test.o
#
mv a.out poisson_1d_multigrid_test
./poisson_1d_multigrid_test > poisson_1d_multigrid_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm poisson_1d_multigrid_test
#
echo "Normal end of execution."
