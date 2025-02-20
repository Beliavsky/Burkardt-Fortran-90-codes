#! /bin/bash
#
gfortran -c -Wall bdf2_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o bdf2_test bdf2_test.o \
  $HOME/lib/bdf2.o \
  $HOME/lib/fsolve.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm bdf2_test.o
#
./bdf2_test > bdf2_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm bdf2_test
#
gnuplot < predator_prey_commands.txt
gnuplot < stiff_commands.txt
#
echo "Normal end of execution."
