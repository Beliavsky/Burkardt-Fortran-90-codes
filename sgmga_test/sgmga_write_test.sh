#! /bin/bash
#
gfortran -c -Wall sgmga_write_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o sgmga_write_test \
  sgmga_write_test.o \
  $HOME/lib/sgmga.o \
  $HOME/lib/sandia_rules.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sgmga_write_test.o
#
./sgmga_write_test > sgmga_write_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sgmga_write_test
#
echo "Normal end of execution."
