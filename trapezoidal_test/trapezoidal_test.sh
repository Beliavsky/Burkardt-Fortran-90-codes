#! /bin/bash
#
gfortran -c -Wall trapezoidal_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o trapezoidal_test trapezoidal_test.o \
  $HOME/lib/trapezoidal.o \
  $HOME/lib/fsolve.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm trapezoidal_test.o
#
./trapezoidal_test > trapezoidal_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm trapezoidal_test
#
gnuplot < predator_prey_commands.txt
gnuplot < stiff_commands.txt
#
echo "Normal end of execution."
