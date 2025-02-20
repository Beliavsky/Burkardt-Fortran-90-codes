#! /bin/bash
#
gfortran -c -Wall slap_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o slap_test slap_test.o $HOME/lib/slap.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm slap_test.o
#
./slap_test > slap_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm slap_test
#
echo "Normal end of execution."
