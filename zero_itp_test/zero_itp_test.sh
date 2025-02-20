#! /bin/bash
#
gfortran -c -Wall zero_itp_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran zero_itp_test.o $HOME/lib/zero_itp.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm zero_itp_test.o
#
mv a.out zero_itp_test
./zero_itp_test > zero_itp_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm zero_itp_test
#
echo "Normal end of execution."
