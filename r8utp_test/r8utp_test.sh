#! /bin/bash
#
gfortran -c -Wall r8utp_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8utp_test r8utp_test.o /$HOME/lib/r8utp.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8utp_test.o
#
./r8utp_test > r8utp_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8utp_test
#
echo "Normal end of execution."
