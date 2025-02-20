#! /bin/bash
#
gfortran -c -Wall sgmga_size_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o sgmga_size_test \
  sgmga_size_test.o \
  $HOME/lib/sgmga.o \
  $HOME/lib/sandia_rules.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sgmga_size_test.o
#
./sgmga_size_test > sgmga_size_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sgmga_size_test
#
echo "Normal end of execution."
