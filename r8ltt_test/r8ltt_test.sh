#! /bin/bash
#
gfortran -c -Wall r8ltt_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8ltt_test r8ltt_test.o $HOME/lib/r8ltt.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8ltt_test.o
#
./r8ltt_test > r8ltt_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8ltt_test
#
echo "Normal end of execution."
