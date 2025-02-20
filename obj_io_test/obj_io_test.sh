#! /bin/bash
#
gfortran -c -Wall obj_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o obj_io_test obj_io_test.o $HOME/lib/obj_io.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm obj_io_test.o
#
./obj_io_test > obj_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm obj_io_test
#
echo "Normal end of execution."
