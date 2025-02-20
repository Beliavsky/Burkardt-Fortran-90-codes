#! /bin/bash
#
gfortran -c -Wall sgmga_vcn_coef_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o sgmga_vcn_coef_test \
  sgmga_vcn_coef_test.o \
  $HOME/lib/sgmga.o \
  $HOME/lib/sandia_rules.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sgmga_vcn_coef_test.o
#
./sgmga_vcn_coef_test > sgmga_vcn_coef_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sgmga_vcn_coef_test
#
echo "Normal end of execution."
