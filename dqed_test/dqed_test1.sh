#! /bin/bash
#
gfortran -c -Wall dqed_test1.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o dqed_test1 dqed_test1.o $HOME/lib/dqed.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm dqed_test1.o
#
./dqed_test1 < dqed_test1_input.txt > dqed_test1.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm dqed_test1
#
echo "Normal end of execution."
