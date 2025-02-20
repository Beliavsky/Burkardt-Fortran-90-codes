#! /bin/bash
#
gfortran -c -Wall r8bto_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o r8bto_test r8bto_test.o $HOME/lib/r8bto.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r8bto_test.o
#
./r8bto_test > r8bto_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r8bto_test
#
echo "Normal end of execution."
