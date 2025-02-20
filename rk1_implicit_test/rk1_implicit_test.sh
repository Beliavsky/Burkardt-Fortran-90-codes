#! /bin/bash
#
gfortran -c -Wall rk1_implicit_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o rk1_implicit_test rk1_implicit_test.o \
  $HOME/lib/rk1_implicit.o \
  $HOME/lib/fsolve.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm rk1_implicit_test.o
#
./rk1_implicit_test > rk1_implicit_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm rk1_implicit_test
#
gnuplot < predator_prey_commands.txt
gnuplot < stiff_commands.txt
#
echo "Normal end of execution."
