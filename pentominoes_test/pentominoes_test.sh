#! /bin/bash
#
gfortran -c -Wall pentominoes_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o pentominoes_test pentominoes_test.o -L$HOME/lib -lpentominoes
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm pentominoes_test.o
#
./pentominoes_test > pentominoes_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm pentominoes_test
#
echo "Normal end of execution."
