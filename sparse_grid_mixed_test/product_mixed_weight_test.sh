#! /bin/bash
#
gfortran -c -Wall product_mixed_weight_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran product_mixed_weight_test.o $HOME/lib/sparse_grid_mixed.o $HOME/lib/sandia_rules.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm product_mixed_weight_test.o
#
mv a.out product_mixed_weight_test
./product_mixed_weight_test > product_mixed_weight_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm product_mixed_weight_test
#
echo "Normal end of execution."
