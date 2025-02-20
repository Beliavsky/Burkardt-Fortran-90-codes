#! /bin/bash
#
gfortran -c -Wall c8poly_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o c8poly_test c8poly_test.o $HOME/lib/c8poly.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm c8poly_test.o
#
./c8poly_test > c8poly_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm c8poly_test
#
echo "Normal end of execution."
