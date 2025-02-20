#! /bin/bash
#
gfortran -c -Wall xy_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran xy_io_test.o -L$HOME/lib -lxy_io
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm xy_io_test.o
#
mv a.out xy_io_test
./xy_io_test > xy_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm xy_io_test
#
echo "Normal end of execution."
