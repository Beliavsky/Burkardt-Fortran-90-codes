#! /bin/bash
#
gfortran -c -Wall rk2_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o rk2_test rk2_test.o $HOME/lib/rk2.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm rk2_test.o
#
./rk2_test > rk2_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm rk2_test
#
gnuplot < predator_prey_rk2_commands.txt
gnuplot < stiff_rk2_commands.txt
#
echo "Normal end of execution."
