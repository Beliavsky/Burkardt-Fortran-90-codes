#! /bin/bash
#
gfortran -c -Wall r8sm_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8sm_test r8sm_test.o $HOME/lib/r8sm.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8sm_test.o
#
./r8sm_test > r8sm_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8sm_test
#
echo "Normal end of execution."
