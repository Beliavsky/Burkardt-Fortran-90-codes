#! /bin/bash
#
gfortran -c -Wall flame_exact_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o flame_exact_test flame_exact_test.o $HOME/lib/flame_exact.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm flame_exact_test.o
#
./flame_exact_test > flame_exact_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm flame_exact_test
#
gnuplot < flame_exact_commands.txt
#
echo "Normal end of execution."
