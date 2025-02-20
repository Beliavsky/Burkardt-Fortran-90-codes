#! /bin/bash
#
gfortran -c -Wall stochastic_heat2d_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o stochastic_heat2d_test stochastic_heat2d_test.o $HOME/lib/stochastic_heat2d.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm stochastic_heat2d_test.o
#
./stochastic_heat2d_test > stochastic_heat2d_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm stochastic_heat2d_test
#
gnuplot < solution_commands.txt
gnuplot < umean_commands.txt
#
echo "Normal end of execution."
