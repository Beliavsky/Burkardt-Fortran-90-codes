#! /bin/bash
#
gfortran -c -Wall steam_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o steam_test steam_test.o $HOME/lib/steam.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm steam_test.o
#
./steam_test > steam_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm steam_test
#
echo "Normal end of execution."
