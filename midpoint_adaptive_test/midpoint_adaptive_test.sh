#! /bin/bash
#
gfortran -c -Wall midpoint_adaptive_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o midpoint_adaptive_test \
  midpoint_adaptive_test.o \
  $HOME/lib/midpoint_adaptive.o \
  $HOME/lib/fsolve.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm midpoint_adaptive_test.o
#
./midpoint_adaptive_test > midpoint_adaptive_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm midpoint_adaptive_test
#
gnuplot < lotka_solution_commands.txt
gnuplot < lotka_phase_commands.txt
gnuplot < lotka_timestep_commands.txt
#
echo "Normal end of execution."
