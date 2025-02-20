#! /bin/bash
#
gfortran -c -Wall comp_next_test.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran comp_next_test.o $HOME/lib/sparse_grid_mixed.o $HOME/lib/sandia_rules.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm comp_next_test.o
#
mv a.out comp_next_test
./comp_next_test > comp_next_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm comp_next_test
#
echo "Normal end of execution."
