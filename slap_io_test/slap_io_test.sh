#! /bin/bash
#
gfortran -c -Wall slap_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o slap_io_test slap_io_test.o $HOME/lib/slap_io.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm slap_io_test.o
#
./slap_io_test > slap_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm slap_io_test
#
echo "Normal end of execution."
