#! /bin/bash
#
gfortran -c -Wall rk2_implicit_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o rk2_implicit_test rk2_implicit_test.o \
  $HOME/lib/rk2_implicit.o \
  $HOME/lib/fsolve.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm rk2_implicit_test.o
#
./rk2_implicit_test > rk2_implicit_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm rk2_implicit_test
#
gnuplot < humps_commands.txt
gnuplot < predator_prey_commands.txt
gnuplot < stiff_commands.txt
#
echo "Normal end of execution."
