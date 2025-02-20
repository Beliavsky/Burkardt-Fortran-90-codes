#! /bin/bash
#
gfortran -c -Wall xyz_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o xyz_io_test xyz_io_test.o $HOME/lib/xyz_io.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm xyz_io_test.o
#
./xyz_io_test > xyz_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm xyz_io_test
#
echo "Normal end of execution."
