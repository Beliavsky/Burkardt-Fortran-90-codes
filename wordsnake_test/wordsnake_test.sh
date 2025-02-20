#! /bin/bash
#
gfortran -c -Wall wordsnake_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o wordsnake_test wordsnake_test.o $HOME/lib/wordsnake.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm wordsnake_test.o
#
./wordsnake_test wordlist.txt > wordsnake_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm wordsnake_test
#
echo "Normal end of execution."
