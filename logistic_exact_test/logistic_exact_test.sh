#! /bin/bash
#
gfortran -c -Wall logistic_exact_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o logistic_exact_test logistic_exact_test.o $HOME/lib/logistic_exact.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm logistic_exact_test.o
#
./logistic_exact_test > logistic_exact_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm logistic_exact_test
#
gnuplot < logistic_exact_commands.txt
#
echo "Normal end of execution."
