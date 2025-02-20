#! /bin/bash
#
gfortran -c -Wall hypersphere_monte_carlo_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran hypersphere_monte_carlo_test.o $HOME/lib/hypersphere_monte_carlo.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm hypersphere_monte_carlo_test.o
#
mv a.out hypersphere_monte_carlo_test
./hypersphere_monte_carlo_test > hypersphere_monte_carlo_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm hypersphere_monte_carlo_test
#
echo "Normal end of execution."
