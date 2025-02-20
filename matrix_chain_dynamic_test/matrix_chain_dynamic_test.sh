#! /bin/bash
#
gfortran -c -Wall matrix_chain_dynamic_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o matrix_chain_dynamic_test matrix_chain_dynamic_test.o $HOME/lib/matrix_chain_dynamic.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm matrix_chain_dynamic_test.o
#
./matrix_chain_dynamic_test > matrix_chain_dynamic_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm matrix_chain_dynamic_test
#
echo "Normal end of execution."
