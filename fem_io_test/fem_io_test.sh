#! /bin/bash
#
gfortran -c -Wall fem_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran fem_io_test.o $HOME/lib/fem_io.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm fem_io_test.o
#
mv a.out fem_io_test
./fem_io_test > fem_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm fem_io_test
#
echo "Normal end of execution."
