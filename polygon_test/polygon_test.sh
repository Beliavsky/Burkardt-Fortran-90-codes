#! /bin/bash
#
gfortran -c -g -Wall polygon_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran polygon_test.o $HOME/lib/polygon.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm polygon_test.o
#
mv a.out polygon_test
./polygon_test > polygon_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm polygon_test
#
echo "Normal end of execution."
