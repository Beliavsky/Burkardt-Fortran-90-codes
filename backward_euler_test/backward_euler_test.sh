#! /bin/bash
#
gfortran -c -Wall backward_euler_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o backward_euler_test backward_euler_test.o \
  $HOME/lib/backward_euler.o \
  $HOME/lib/fsolve.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm backward_euler_test.o
#
./backward_euler_test > backward_euler_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm backward_euler_test
#
gnuplot < predator_prey_commands.txt
gnuplot < stiff_commands.txt
#
echo "Normal end of execution."
