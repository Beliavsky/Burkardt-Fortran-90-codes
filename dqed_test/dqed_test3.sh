#! /bin/bash
#
gfortran -c -Wall dqed_test3.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o dqed_test3 dqed_test3.o $HOME/lib/dqed.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm dqed_test3.o
#
./dqed_test3 > dqed_test3.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm dqed_test3
#
echo "Normal end of execution."
