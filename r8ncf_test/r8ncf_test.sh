#! /bin/bash
#
gfortran -c -Wall r8ncf_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8ncf_test r8ncf_test.o $HOME/lib/r8ncf.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8ncf_test.o
#
./r8ncf_test > r8ncf_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8ncf_test
#
echo "Normal end of execution."
