#! /bin/bash
#
gfortran -c -Wall rk1_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o rk1_test rk1_test.o $HOME/lib/rk1.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm rk1_test.o
#
./rk1_test > rk1_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm rk1_test
#
echo "Generate png image."
#
gnuplot < rk1_humps_commands.txt
#
echo "Normal end of execution."
