#! /bin/bash
#
gfortran -c -Wall sgmga_index_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o sgmga_index_test \
  sgmga_index_test.o \
  $HOME/lib/sgmga.o \
  $HOME/lib/sandia_rules.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sgmga_index_test.o
#
./sgmga_index_test > sgmga_index_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sgmga_index_test
#
echo "Normal end of execution."
