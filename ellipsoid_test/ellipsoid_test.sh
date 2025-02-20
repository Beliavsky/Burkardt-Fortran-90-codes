#! /bin/bash
#
gfortran -c -Wall ellipsoid_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o ellipsoid_test ellipsoid_test.o $HOME/lib/ellipsoid.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm ellipsoid_test.o
#
./ellipsoid_test > ellipsoid_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm ellipsoid_test
#
echo "Normal end of execution."
