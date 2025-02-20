#! /bin/bash
#
gfortran -c -Wall sigmoid_derivative_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o sigmoid_derivative_test sigmoid_derivative_test.o $HOME/lib/sigmoid_derivative.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sigmoid_derivative_test.o
#
./sigmoid_derivative_test > sigmoid_derivative_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sigmoid_derivative_test
#
gnuplot sigmoid_derivative_0_commands.txt
gnuplot sigmoid_derivative_1_commands.txt
gnuplot sigmoid_derivative_2_commands.txt
gnuplot sigmoid_derivative_3_commands.txt
#
echo "Normal end of execution."
