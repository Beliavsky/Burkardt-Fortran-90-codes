#! /bin/bash
#
gfortran -c -Wall dqed_test2.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o dqed_test2 dqed_test2.o $HOME/lib/dqed.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm dqed_test2.o
#
./dqed_test2 > dqed_test2.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm dqed_test2
#
echo "Normal end of execution."
