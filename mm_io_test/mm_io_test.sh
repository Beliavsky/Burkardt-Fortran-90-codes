#! /bin/bash
#
gfortran -c -Wall mm_io_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o mm_io_test mm_io_test.o $HOME/lib/mm_io.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm mm_io_test.o
#
./mm_io_test > mm_io_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm mm_io_test
#
echo "Normal end of execution."
