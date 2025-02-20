#! /bin/bash
#
gfortran -c -Wall biharmonic_exact_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o biharmonic_exact_test \
  biharmonic_exact_test.o \
  $HOME/lib/biharmonic_exact.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm biharmonic_exact_test.o
#
./biharmonic_exact_test > biharmonic_exact_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm biharmonic_exact_test
#
echo "Normal end of execution."
