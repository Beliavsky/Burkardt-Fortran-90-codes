#! /bin/bash
#
gfortran -c -Wall r8ccs_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8ccs_test r8ccs_test.o $HOME/lib/r8ccs.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8ccs_test.o
#
./r8ccs_test > r8ccs_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8ccs_test
#
echo "Normal end of execution."
