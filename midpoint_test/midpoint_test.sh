#! /bin/bash
#
gfortran -c -Wall midpoint_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o midpoint_test midpoint_test.o \
  $HOME/lib/midpoint.o \
  $HOME/lib/fsolve.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm midpoint_test.o
#
./midpoint_test > midpoint_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm midpoint_test
#
gnuplot < humps_commands.txt
gnuplot < predator_prey_commands.txt
gnuplot < stiff_commands.txt
#
echo "Normal end of execution."
