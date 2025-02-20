#! /bin/bash
#
gfortran -c -Wall porous_medium_exact_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o porous_medium_exact_test porous_medium_exact_test.o \
  $HOME/lib/porous_medium_exact.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm porous_medium_exact_test.o
#
./porous_medium_exact_test > porous_medium_exact_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm porous_medium_exact_test
#
echo "Normal end of execution."
