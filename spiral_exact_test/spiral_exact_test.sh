#! /bin/bash
#
gfortran -c -Wall spiral_exact_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o spiral_exact_test spiral_exact_test.o $HOME/lib/spiral_exact.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm spiral_exact_test.o
#
./spiral_exact_test > spiral_exact_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm spiral_exact_test
#
gnuplot < spiral_commands.txt
#
echo "Normal end of execution."
