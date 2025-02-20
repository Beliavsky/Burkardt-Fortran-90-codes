#! /bin/bash
#
gfortran -c -Wall slatec_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o slatec_test slatec_test.o $HOME/lib/slatec.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm slatec_test.o
#
./slatec_test > slatec_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm slatec_test
#
echo "Normal end of execution."
