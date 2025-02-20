#! /bin/bash
#
gfortran -c -Wall line_monte_carlo_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran line_monte_carlo_test.o $HOME/lib/line_monte_carlo.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm line_monte_carlo_test.o
#
mv a.out line_monte_carlo_test
./line_monte_carlo_test > line_monte_carlo_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm line_monte_carlo_test
#
echo "Normal end of execution."
