#! /bin/bash
#
gfortran -c -Wall sgmga_size_table.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran -o sgmga_size_table \
  sgmga_size_table.o \
  $HOME/lib/sgmga.o \
  $HOME/lib/sandia_rules.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm sgmga_size_table.o
#
./sgmga_size_table > sgmga_size_table.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm sgmga_size_table
#
echo "Normal end of execution."
