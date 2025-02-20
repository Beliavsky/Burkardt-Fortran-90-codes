#! /bin/bash
#
gfortran -c -Wall r8sto_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8sto_test r8sto_test.o $HOME/lib/r8sto.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8sto_test.o
#
./r8sto_test > r8sto_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8sto_test
#
echo "Normal end of execution."
