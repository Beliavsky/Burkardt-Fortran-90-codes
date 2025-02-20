#! /bin/bash
#
gfortran -c -Wall r8pbu_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8pbu_test r8pbu_test.o $HOME/lib/r8pbu.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8pbu_test.o
#
./r8pbu_test > r8pbu_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8pbu_test
#
echo "Normal end of execution."
