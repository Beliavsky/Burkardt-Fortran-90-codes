#! /bin/bash
#
gfortran -c -Wall dqed_test5.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o dqed_test5 dqed_test5.o $HOME/lib/dqed.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm dqed_test5.o
#
./dqed_test5 > dqed_test5.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm dqed_test5
#
echo "Normal end of execution."
