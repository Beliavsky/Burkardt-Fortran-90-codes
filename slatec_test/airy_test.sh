#! /bin/bash
#
gfortran -c -Wall airy_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o airy_test airy_test.o $HOME/lib/slatec.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm airy_test.o
#
./airy_test > airy_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm airy_test
#
echo "Normal end of execution."
