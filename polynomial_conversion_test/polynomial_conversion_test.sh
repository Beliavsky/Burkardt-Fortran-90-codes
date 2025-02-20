#! /bin/bash
#
gfortran -c -cpp -Wall polynomial_conversion_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o polynomial_conversion_test polynomial_conversion_test.o $HOME/lib/polynomial_conversion.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm polynomial_conversion_test.o
#
./polynomial_conversion_test > polynomial_conversion_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm polynomial_conversion_test
#
echo "Normal end of execution."
