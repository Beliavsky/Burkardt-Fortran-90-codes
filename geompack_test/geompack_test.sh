#! /bin/bash
#
gfortran -c -Wall geompack_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o geompack_test geompack_test.o $HOME/lib/geompack.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm geompack_test.o
#
./geompack_test > geompack_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm geompack_test
#
echo "Normal end of execution."
