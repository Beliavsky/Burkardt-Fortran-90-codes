#! /bin/bash
#
gfortran -c -Wall stochastic_rk_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o stochastic_rk_test stochastic_rk_test.o $HOME/lib/stochastic_rk.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm stochastic_rk_test.o
#
./stochastic_rk_test > stochastic_rk_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm stochastic_rk_test
#
echo "Normal end of execution."
