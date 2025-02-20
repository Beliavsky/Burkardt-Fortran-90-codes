#! /bin/bash
#
gfortran -c -Wall fisher_exact_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o fisher_exact_test fisher_exact_test.o $HOME/lib/fisher_exact.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm fisher_exact_test.o
#
./fisher_exact_test > fisher_exact_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm fisher_exact_test
#
echo "Normal end of execution."
