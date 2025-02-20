#! /bin/bash
#
gfortran -c -g -Wall triangle_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o triangle_test triangle_test.o \
  $HOME/lib/triangle.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm triangle_test.o
#
./triangle_test > triangle_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm triangle_test
#
echo "Normal end of execution."
